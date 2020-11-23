#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jpeglib.h"
#include <setjmp.h>
#include <math.h>
#include "../simple-jpeg/importexport.c"
#include <mpi.h>
#include <sys/time.h>

typedef struct
{
    float** image_data;
    int m; // height
    int n; // width
}
image;

void allocate_image(image *u, int m, int n);
void deallocate_image(image *u);
void convert_jpeg_to_image(const unsigned char* image_chars, image *u);
void convert_image_to_jpeg(const image *u, unsigned char* image_chars);
void iso_diffusion_denoising_parallel(image *u, image *u_bar, float kappa, int iters, int my_n_rank, int my_m_rank, int num_n_procs, int num_m_procs, int my_rank);

int get_optimal_m_axis_process_count(int workers, int m, int y);
double get_time();
double get_time_passed(double previous);

int main(int argc, char *argv[])
{
    // master variables
    int *external_m_starts, *external_n_starts, *external_ms, *external_ns;
    unsigned char *image_chars, *transfer_chars;
    char *input_jpeg_filename, *output_jpeg_filename;
    double my_time;

    // universal variables
    int m, n, c, iters, my_m, my_n, my_rank, num_procs;
    int my_n_rank, my_m_rank, num_n_procs, num_m_procs;
    float kappa;
    image u, u_bar;
    MPI_Status my_stat;

    // worker variables
    unsigned char *my_image_chars;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

    if(my_rank == 0)
    {
        if(argc < 5)
        {
            printf("%u/5 arguments provided, running with default parameters\n", argc);
            kappa = 0.2f;
            iters = 100;
            input_jpeg_filename = "mona_lisa_noisy.jpg";
            output_jpeg_filename = "mona_lisa_quiet.jpg";
        }
        else
        {
            kappa = atof(argv[1]);
            iters = (int)argv[2];
            input_jpeg_filename = argv[3];
            output_jpeg_filename = argv[4];
        }
        import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
        printf("Image properties:  width - %u, height - %u, c - %u\n", n, m, c);
        printf("Preparing to denoise\n");
        printf("********************************************\n");
        fflush(stdout);
    }

    // we give the worker threads enough info that they can figure out the rest themselves
    MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&kappa, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&iters, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&num_n_procs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // it is better that the threads find these numbers individually, so as to save on bandwidth
    num_m_procs = get_optimal_m_axis_process_count(num_procs, m, n);
    num_n_procs = num_procs / num_m_procs;
    my_n_rank = my_rank % num_n_procs;
    my_m_rank = my_rank / num_n_procs;
    my_n = n / num_n_procs;
    if(my_n_rank < n % num_n_procs)
        my_n++;
    my_m = m / num_m_procs;
    if(my_m_rank < m % num_m_procs)
        my_m++;

    printf("Thread %u has indices n:%u m:%u, %ux%u pixels\n", my_rank, my_n_rank, my_m_rank, my_n, my_m);
    fflush(stdout);

    // 0 finds and allocates the parameters it needs to split the jpg
    if(my_rank == 0)
    {
        int export_n_rank, export_m_rank;
        external_ms = (int*) malloc(num_procs * sizeof(int));
        external_ns = (int*) malloc(num_procs * sizeof(int));
        external_m_starts = (int*) malloc(num_procs * sizeof(int));
        external_n_starts = (int*) malloc(num_procs * sizeof(int));
        transfer_chars = (unsigned char*) malloc(my_n * my_m * sizeof(unsigned char)); // dimensions are ok, because 0 has the largest area

        for(int i = 1; i < num_procs; i++)
        {
             export_n_rank = i % num_n_procs;
             export_m_rank = i / num_n_procs;
             external_ns[i] = n / num_n_procs;
             if(export_n_rank < n % num_n_procs)
                 external_ns[i]++;
             external_ms[i] = m / num_m_procs;
             if(export_m_rank < m % num_m_procs)
                 external_ms[i]++;
             external_n_starts[i] = export_n_rank * (n / num_n_procs) + fmin(export_n_rank, n % num_n_procs);
             external_m_starts[i] = export_m_rank * (m / num_m_procs) + fmin(export_m_rank, m % num_m_procs);
        }
    }

    for(int i = 0; i < c; i++)
    {
        // 0 distributes the parts evenly
        if(my_rank == 0)
        {
            for(int ii = 1; ii < num_procs; ii++)
            {
                 for(int iii = 0; iii < external_ms[ii]; iii++)
                     for(int iv = 0; iv < external_ns[ii]; iv++)
                         transfer_chars[iii * external_ns[ii] + iv] = image_chars[((iii + external_m_starts[ii]) * n + external_n_starts[ii] + iv) * c + i];

                 MPI_Send(transfer_chars, external_ms[ii] * external_ns[ii], MPI_UNSIGNED_CHAR, ii, 4, MPI_COMM_WORLD);
            }

            // 0 then takes its own part of the jpg
            my_image_chars = (unsigned char*) malloc(my_n * my_m * sizeof(unsigned char));
            for(int i = 0; i < my_m; i++)
                for(int ii = 0; ii < my_n; ii++)
                    my_image_chars[i * my_n + ii] = image_chars[(i * n + ii) * c + i];
            printf("********************************************\n");
        }
        else    // everyone else waits to receive their share of the jpg
        {
            my_image_chars = (unsigned char*) malloc(my_n * my_m * sizeof(unsigned char));
            MPI_Recv(my_image_chars, my_m * my_n, MPI_UNSIGNED_CHAR, 0, 4, MPI_COMM_WORLD, &my_stat);
        }

        // all processes convert and denoise their piece of the jpg
        allocate_image (&u, my_m, my_n);
        allocate_image (&u_bar, my_m, my_n);
        convert_jpeg_to_image (my_image_chars, &u);
        if(my_rank == 0)
        {
            printf("Denoising image, component %u\n", i);
            fflush(stdout);
            my_time = get_time();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        iso_diffusion_denoising_parallel (&u, &u_bar, kappa, iters, my_n_rank, my_m_rank, num_n_procs, num_m_procs, my_rank);
        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank == 0)
        {
            my_time = get_time_passed(my_time);
            printf("Image denoised, %fs passed\n", my_time);
            printf("********************************************\n");
            fflush(stdout);
        }
        convert_image_to_jpeg(&u_bar, my_image_chars);

        // 0 recieves the pieces back from the other processes
        if(my_rank == 0)
        {
            // inserts own data first
            for(int ii = 0; ii < my_m; ii++)
                for(int iii = 0; iii < my_n; iii++)
                    image_chars[(ii * n + iii) * c + i] = my_image_chars[ii * my_n + iii];

            // then the rest
            for(int ii = 1; ii < num_procs; ii++)
            {
                MPI_Recv(transfer_chars, external_ms[ii] * external_ns[ii], MPI_UNSIGNED_CHAR, ii, 0, MPI_COMM_WORLD, &my_stat);

                for(int iii = 0; iii < external_ms[ii]; iii++)
                    for(int iv = 0; iv < external_ns[ii]; iv++)
                        image_chars[((iii + external_m_starts[ii]) * n + external_n_starts[ii] + iv) * c + i] = transfer_chars[iii * external_ns[ii] + iv];

            }
        }
        else // everyone else waits to send their share of the jpg
        {
            MPI_Send(my_image_chars, my_m * my_n, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
        }
    }

    // finished; export, deallocate and return
    if(my_rank==0)
    {
        export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
        free(external_ms);
        free(external_ns);
        free(external_m_starts);
        free(external_n_starts);
        free(image_chars);
        free(transfer_chars);
        printf("Denoise complete!\n");
    }
    else
    {
        free(my_image_chars);
    }
    deallocate_image (&u);
    deallocate_image (&u_bar);
    MPI_Finalize ();
    return 0;
}

void iso_diffusion_denoising_parallel(image *u, image *u_bar, float kappa, int iters, int my_n_rank, int my_m_rank, int num_n_procs, int num_m_procs, int my_rank)
{    
    int m = u->m, n = u->n;
    float *from_beyond_min_n, *from_beyond_min_m, *from_beyond_max_n, *from_beyond_max_m;
    float *to_beyond_min_n, *to_beyond_min_m, *to_beyond_max_n, *to_beyond_max_m;
    MPI_Status my_stat;
    MPI_Request my_reqtm, my_reqbm, my_reqtn, my_reqbn;

    // allocating space for external data
    if(my_n_rank > 0)
    {
        from_beyond_min_n = (float*) malloc(m * sizeof(float));
        to_beyond_min_n = (float*) malloc(m * sizeof(float));
    }
    if(my_m_rank > 0)
    {
        from_beyond_min_m = (float*) malloc(n * sizeof(float));
        to_beyond_min_m = (float*) malloc(n * sizeof(float));
    }
    if(my_n_rank < num_n_procs - 1)
    {
        from_beyond_max_n = (float*) malloc(m * sizeof(float));
        to_beyond_max_n = (float*) malloc(m * sizeof(float));
    }
    if(my_m_rank < num_m_procs - 1)
    {
        from_beyond_max_m = (float*) malloc(n * sizeof(float));
        to_beyond_max_m = (float*) malloc(n * sizeof(float));
    }

    for(int i = 0; i < iters; i++)
    {
        // initialize u barred
        for(int ii = 0; ii < m; ii++)
            for(int iii = 0; iii < n; iii++)
                u_bar->image_data[ii][iii] = 0.0f;

        // export, import and add external edges
        if(my_n_rank > 0)
        {
            for(int ii = 0; ii < m; ii++)
                to_beyond_min_n[ii] = u->image_data[ii][0];
            MPI_Isend(to_beyond_min_n, m, MPI_FLOAT, my_rank - 1, 0, MPI_COMM_WORLD, &my_reqbn);
            MPI_Recv(from_beyond_min_n, m, MPI_FLOAT, my_rank - 1, 1, MPI_COMM_WORLD, &my_stat);
            for(int ii = 0; ii < m; ii++)
                u_bar->image_data[ii][0] += from_beyond_min_n[ii] - u->image_data[ii][0];
        }
        if(my_n_rank < num_n_procs - 1)
        {
            for(int ii = 0; ii < m; ii++)
                to_beyond_max_n[ii] = u->image_data[ii][n - 1];
            MPI_Isend(to_beyond_max_n, m, MPI_FLOAT, my_rank + 1, 1, MPI_COMM_WORLD, &my_reqtn);
            MPI_Recv(from_beyond_max_n, m, MPI_FLOAT, my_rank + 1, 0, MPI_COMM_WORLD, &my_stat);
            for(int ii = 0; ii < m; ii++)
                u_bar->image_data[ii][n - 1] += from_beyond_max_n[ii] - u->image_data[ii][n - 1];
        }
        if(my_m_rank > 0)
        {
            for(int ii = 0; ii < n; ii++)
                to_beyond_min_m[ii] = u->image_data[0][ii];
            MPI_Isend(to_beyond_min_m, n, MPI_FLOAT, my_rank - num_n_procs, 2, MPI_COMM_WORLD, &my_reqbm);
            MPI_Recv(from_beyond_min_m, n, MPI_FLOAT, my_rank - num_n_procs, 3, MPI_COMM_WORLD, &my_stat);
            for(int ii = 0; ii < n; ii++)
                u_bar->image_data[0][ii] += from_beyond_min_m[ii] - u->image_data[0][ii];
        }
        if(my_m_rank < num_m_procs - 1)
        {
            for(int ii = 0; ii < n; ii++)
                to_beyond_max_m[ii] = u->image_data[m - 1][ii];
            MPI_Isend(to_beyond_max_m, n, MPI_FLOAT, my_rank + num_n_procs, 3, MPI_COMM_WORLD, &my_reqtm);
            MPI_Recv(from_beyond_max_m, n, MPI_FLOAT, my_rank + num_n_procs, 2, MPI_COMM_WORLD, &my_stat);
            for(int ii = 0; ii < n; ii++)
                u_bar->image_data[m - 1][ii] += from_beyond_max_m[ii] - u->image_data[m - 1][ii];
        }

        // -1 m
        for(int ii = 1; ii < m; ii++)
            for(int iii = 0; iii < n; iii++)
                u_bar->image_data[ii][iii] += u->image_data[ii - 1][iii] - u->image_data[ii][iii];
        // -1 n
        for(int ii = 0; ii < m; ii++)
            for(int iii = 1; iii < n; iii++)
                u_bar->image_data[ii][iii] += u->image_data[ii][iii - 1] - u->image_data[ii][iii];
        // +1 m
        for(int ii = 0; ii < m - 1; ii++)
            for(int iii = 0; iii < n; iii++)
                u_bar->image_data[ii][iii] += u->image_data[ii + 1][iii] - u->image_data[ii][iii];
        // +1 n
        for(int ii = 0; ii < m; ii++)
            for(int iii = 0; iii < n - 1; iii++)
                u_bar->image_data[ii][iii] += u->image_data[ii][iii + 1] - u->image_data[ii][iii];
        // finishing equation
        for(int ii = 0; ii < m; ii++)
            for(int iii = 0; iii < n; iii++)
                u_bar->image_data[ii][iii] = u_bar->image_data[ii][iii] * kappa + u->image_data[ii][iii];

        // not last iteration
        if(i < iters - 1)
        {
            // move u bar to u
            for(int ii = 0; ii < m; ii++)
                for(int iii = 0; iii < n; iii++)
                    u->image_data[ii][iii] = u_bar->image_data[ii][iii];

            // make sure the send buffer is ready before starting a new iteration
            if(my_n_rank > 0)
                MPI_Wait(&my_reqbn, &my_stat);
            if(my_m_rank > 0)
                MPI_Wait(&my_reqbm, &my_stat);
            if(my_n_rank < num_n_procs - 1)
                MPI_Wait(&my_reqtn, &my_stat);
            if(my_m_rank < num_m_procs - 1)
                MPI_Wait(&my_reqtm, &my_stat);
        }
        else
        {
            // last iteration, deallocate and return
            if(my_n_rank > 0)
            {
                free(from_beyond_min_n);
                free(to_beyond_min_n);
            }
            if(my_m_rank > 0)
            {
                free(from_beyond_min_m);
                free(to_beyond_min_m);
            }
            if(my_n_rank < num_n_procs - 1)
            {
                free(from_beyond_max_n);
                free(to_beyond_max_n);
            }
            if(my_m_rank < num_m_procs - 1)
            {
                free(from_beyond_max_m);
                free(to_beyond_max_m);
            }
        }
    }
}

void allocate_image(image *u, int m, int n)
{
    u->image_data = (float**) malloc(m * sizeof(float*));
    for(int i = 0; i < m; i++)
    {
        u->image_data[i] = (float*) malloc(n * sizeof(float));
    }

    u->m = m;
    u->n = n;
}

void deallocate_image(image *u)
{
    for(int i = 0; i < u->m; i++)
    {
        free(u->image_data[i]);
    }
    free(u->image_data);
}

void convert_jpeg_to_image(const unsigned char *image_chars, image *u)
{
    double mult = 1.0 / 255.0;
    for(int i = 0; i < u->m; i++)
    {
        for(int ii = 0; ii < u->n; ii++)
        {
            u->image_data[i][ii] = (double) image_chars[i * u->n + ii] * mult;
        }
    }
}

void convert_image_to_jpeg(const image *u, unsigned char *image_chars)
{
    for(int i = 0; i < u->m; i++)
    {
        for(int ii = 0; ii < u->n; ii++)
        {
            image_chars[i * u->n + ii] = (char)(u->image_data[i][ii] * 255);
        }
    }
}

int get_optimal_m_axis_process_count(int workers, int m, int y)
{
    int x = 1, diff = workers;
    for(int i = 1; i < workers; i++)
    {
        if(workers % i == 0 && abs((workers / i) - i) < diff)
        {
            diff = abs((workers / i) - i);
            x = i;
        }
    }
    if(m > y)
        return workers / x;
    return x;
}

double get_time()
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 0.000000001l;
}

double get_time_passed(double previous)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 0.000000001l - previous;
}
