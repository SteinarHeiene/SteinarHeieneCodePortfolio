#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jpeglib.h"

#include <setjmp.h>
#include "../simple-jpeg/importexport.c"

typedef struct
{
float** image_data;
int m; /* # pixels in vertical-direction */
int n; /* # pixels in horizontal-direction */
}
image;

void allocate_image(image *u, int m, int n);
void deallocate_image(image *u);
void convert_jpeg_to_image(const unsigned char* image_chars, image *u);
void convert_image_to_jpeg(const image *u, unsigned char* image_chars);
void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters);

int main(int argc, char *argv[])
{
    int m, n, c, iters;
    float kappa;
    image u, u_bar;
    unsigned char *image_chars;
    char *input_jpeg_filename, *output_jpeg_filename;

    if(argc != 4)
    {
        printf("4 arguments not provided, running with default paramaters\n");
        kappa = 0.2f;
        iters = 100;
        input_jpeg_filename = "mona_lisa_noisy.jpg";
        output_jpeg_filename = "mona_lisa_quiet.jpg";
    }
    else
    {
        kappa = atof(argv[0]);
        iters = (int)argv[1];
        input_jpeg_filename = argv[2];
        output_jpeg_filename = argv[3];
    }

    printf("Reading file: %s ", input_jpeg_filename);
    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    printf("COMPLETE\nWidth: %u\nHeight: %u\nComponents: %u\n", n, m, c);

    allocate_image(&u, m, n);
    allocate_image(&u_bar, m, n);

    convert_jpeg_to_image(image_chars, &u);
    printf("\nDenoising image\nKappa: %f\nIterations: %u\n", kappa, iters);
    iso_diffusion_denoising(&u, &u_bar, kappa, iters);
    printf("COMPLETE\n");

    convert_image_to_jpeg(&u_bar, image_chars);
    export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);

    deallocate_image (&u);
    deallocate_image (&u_bar);

    printf("Application terminating successfully\n");

    return 0;
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

void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters)
{
    int m = u->m, n = u->n;
    for(int i = 0; i < iters; i++)
    {
        // initialize u barred
        for(int ii = 0; ii < m; ii++)
            for(int iii = 0; iii < n; iii++)
                u_bar->image_data[ii][iii] = 0.0f;
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

        if(i < iters - 1)
        {
            // moving u barred to u. We skip edges since they remain unchanged
            for(int ii = 1; ii < m - 1; ii++)
            {
                for(int iii = 1; iii < n - 1; iii++)
                {
                    u->image_data[ii][iii] = u_bar->image_data[ii][iii];
                }
            }
        }
    }
}
