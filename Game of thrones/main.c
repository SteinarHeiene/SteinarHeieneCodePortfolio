#include <stdio.h>
#include <stdlib.h>

void readFile(char* fileName);
void runSimulations();
int permutate(int *collection, int *override);
void printResults();
void cleanUp();

char eNames[32][16];
char charNames[28][16] =
{
    "Jon", "Daenarys", "Tyrion", "Arya", "Cersei", "Sansa", "Bran", "Varys", "Jaime", "Theon", "Yara",
    "Euron", "Samwell", "Gilly", "Brienne", "Davos", "Qyburn", "Grey Worm", "Missandei", "Jorah",
    "The Hound", "Bronn", "Melisandre", "Tormund", "Beric", "Gendry", "Meera", "Night king"
};


int **entries; // alive 0, dead 1, irellevant alive 2, irellevant dead 3
int *solution, *eFirsts, *eLasts, *eThrones;
int eCount, cFirst, cLast;

//statistics:
int *timesWon, *timesDrawn;
int **timesDeadIfEWon, **timesLivedIfEWon, **timesThronedIfEWon, **timesLastIfEWon;

int main()
{   
    readFile("readfile.txt");
    timesWon = (int*) calloc(eCount, sizeof(int));
    timesDrawn = (int*) calloc(eCount, sizeof(int));

    timesDeadIfEWon = (int**) malloc(eCount * sizeof(int*));
    timesLivedIfEWon = (int**) malloc(eCount * sizeof(int*));
    timesThronedIfEWon = (int**) malloc(eCount * sizeof(int*));
    timesLastIfEWon = (int**) malloc(eCount * sizeof(int*));
    for(int i = 0; i < eCount; i++)
    {
        timesDeadIfEWon[i] = (int*) calloc(28, sizeof(int));
        timesLivedIfEWon[i] = (int*) calloc(28, sizeof(int));
        timesThronedIfEWon[i] = (int*) calloc(28, sizeof(int));
        timesLastIfEWon[i] = (int*) calloc(28, sizeof(int));
    }

    runSimulations();
    writeResults();
    //printResults();
    cleanUp();

    return 0;
}

void runSimulations()
{
    printf("Setting up simulation environment...");
    // base, simulation, lasts, thrones
    int *simulation, *fPoints, *fbPoints, *fblPoints, *fbltPoints;
    int sLast, expectedLoops = 1, loop = 0, simCount = 0;

    simulation = (int*) malloc(28 * sizeof(int));
    fPoints = (int*) malloc(eCount * sizeof(int));
    fbPoints = (int*) malloc(eCount * sizeof(int));
    fblPoints = (int*) malloc(eCount * sizeof(int));
    fbltPoints = (int*) malloc(eCount * sizeof(int));

    for(int i = 0; i < 28; i++)
    {
        switch (solution[i])
        {
        case 0:
            expectedLoops *= 2;
        case 1:
            simulation[i] = solution[i];
        break;
        }
    }

    // tally determined firsts points
    for(int i = 0; i < eCount; i++)
    {
        if(eFirsts[i] == cFirst)
        {
            fPoints[i] = 3;
        }
        else
        {
            fPoints[i] = 0;
        }
    }

    printf("Complete\nRunning simulations.........");
    // for each permutation of deaths
    do
    {
        if(loop % 5000 == 0)
        {
            float percent = ((float)loop / (float)expectedLoops) * 100;

            if(percent >= 10.0f)
            {
                printf("\b");
            }
            if(percent >= 100.0f)
            {
                printf("\b");
            }

            printf("\b\b\b\b\b%.2f%%", percent);
        }
        loop++;

        // tally base simulation points
        for(int i = 0; i < eCount; i++)
        {
            fbPoints[i] = fPoints[i];
            for(int u = 0; u < 28; u++)
            {
                if(simulation[u] == entries[i][u])
                {
                    fbPoints[i]++;
                }
            }
        }

        // go through every iteration of lasts
        for(int i = 0; i < 28; i++)
        {
            // as long as the character is dead in the simulation but not in the show
            // if noone died in the simulation use the show last
            if(simulation[i] == 1 && solution[i] == 0 || (loop == 1 && i == sLast))
            {
                // tally last points
                for(int u = 0; u < eCount; u++)
                {
                    fblPoints[u] = fbPoints[u];
                    if(eLasts[u] == i)
                    {
                        fblPoints[u] += 3;
                    }
                }

                // go through every iteration of thrones
                for(int u = 0; u < 28; u++)
                {
                    // as long as the character is alive in the simulation
                    if(simulation[u] == 0)
                    {
                        int highScore = 0;
                        int winnerCount = 0;
                        // tally thrones points and find winners
                        for(int y = 0; y < eCount; y++)
                        {
                            fbltPoints[y] = fblPoints[y];
                            if(eThrones[y] == u)
                            {
                                fbltPoints[y] += 3;
                            }

                            if(fbltPoints[y] > highScore)
                            {
                                highScore = fbltPoints[y];
                                winnerCount = 1;
                            }
                            else if(fbltPoints[y] == highScore)
                            {
                                winnerCount++;
                            }
                        }

                        // update statistics
                        if(winnerCount > 1)
                        {
                            for(int y = 0; y < eCount; y++)
                            {
                                if(fbltPoints[y] == highScore)
                                {
                                    timesDrawn[y]++;
                                }
                            }
                        }
                        else
                        {
                            for(int y = 0; y < eCount; y++)
                            {
                                if(fbltPoints[y] == highScore)
                                {
                                    timesWon[y]++;
                                    for(int t = 0; t < 28; t++)
                                    {
                                        if(simulation[t] == 0)
                                            timesLivedIfEWon[y][t]++;
                                        else
                                            timesDeadIfEWon[y][t]++;
                                    }
                                    timesLastIfEWon[y][i]++;
                                    timesThronedIfEWon[y][u]++;
                                }
                            }
                        }

                        simCount++;
                    }
                }
            }
        }
    }
    while(permutate(simulation, solution) == 0);

    printf("\b\b\b\b\b\bComplete: %u simulations\n\n", simCount);

    free(simulation);
    free(fPoints);
    free(fbPoints);
    free(fbPoints);
    free(fblPoints);
    free(fbltPoints);
}

void readFile(char *fileName)
{
    FILE *readFile = fopen(fileName, "r");
    int readInt;
    char readString[32];

    printf("Reading file %s...", fileName);
    if(readFile == NULL)
    {
        printf("Could not read file: %s\n", fileName);
        exit(1);
    }

    fscanf(readFile, "%u", &cFirst);
    fscanf(readFile, "%u", &cLast);
    fscanf(readFile, "%u", &eCount);
    solution = (int*) malloc(28 * sizeof(int));
    for(int i = 0; i < 28; i++)
    {
        fscanf(readFile, "%u", &solution[i]);
    }

    entries = (int**) malloc(eCount * sizeof(int*));
    eFirsts = (int*) malloc(eCount * sizeof(int));
    eLasts = (int*) malloc(eCount * sizeof(int));
    eThrones = (int*) malloc(eCount * sizeof(int));

    for(int i = 0; i < eCount; i++)
    {
        entries[i] = (int*) malloc(28 * sizeof(int));
        fscanf(readFile, "%s", &eNames[i]);

        for(int u = 0; u < 28; u++)
        {
            fscanf(readFile, "%u", &entries[i][u]);
        }
        fscanf(readFile, "%u", &eFirsts[i]);
        fscanf(readFile, "%u", &eLasts[i]);
        fscanf(readFile, "%u", &eThrones[i]);
    }
    printf("Complete\n\n*****Entries:*****\n");
    for(int i = 0; i < eCount; i++)
    {
        printf("%s: ", eNames[i]);
        for(int u = 0; u < 28; u++)
        {
            printf("%u", entries[i][u]);
        }
        printf(" %u %u %u\n", eFirsts[i], eLasts[i], eThrones[i]);
    }
    printf("\n");
}

int permutate(int* collection, int *override)
{
    for(int i = 0; i < 28; i++)
    {
        if(collection[i] == 0)
        {
            collection[i] = 1;
            return 0;
        }
        else if(override[i] == 0)
        {
            collection[i] = 0;
        }
    }
    return 1;
}

void writeResults()
{
    FILE *outfile;
    int c;
    outfile = fopen("writefile.txt", "w+");
    for(int i = 0; i < eCount; i++)
    {
        fprintf(outfile, "%s\n", eNames[i]);

        for(int u = 0; u < 28; u++)
        {
            fprintf(outfile, "      %u\n", timesDeadIfEWon[i][u]);
            fprintf(outfile, "      %u\n", timesLastIfEWon[i][u]);
            fprintf(outfile, "      %u\n", timesThronedIfEWon[i][u]);
        }
        fprintf(outfile, "   %u\n", timesWon[i]);
        fprintf(outfile, "   %u\n", timesDrawn[i]);
    }
    fclose(outfile);
}

void printResults()
{
    printf("***RESULTS:***\n");
    for(int i = 0; i < eCount; i++)
    {
        printf("%s:\n", eNames[i]);
        printf("Times won: %u Times tied: %u\n", timesWon[i], timesDrawn[i]);
        printf("Character specifics:\n");
        for(int u = 0; u < 28; u++)
        {
            printf("  %s:", charNames[u]);
            printf(" died: %u", timesDeadIfEWon[i][u]);
            printf(" died last: %u", timesLastIfEWon[i][u]);
            printf(" lived: %u", timesLivedIfEWon[i][u]);
            printf(" throned: %u\n", timesThronedIfEWon[i][u]);
        }
        fflush(stdout);
        getchar();
    }
}

void cleanUp()
{
    for(int i = 0; i < eCount; i++)
    {
        free(entries[i]);
        free(timesDeadIfEWon[i]);
        free(timesLivedIfEWon[i]);
        free(timesThronedIfEWon[i]);
        free(timesLastIfEWon[i]);
    }
    free(entries);
    free(eFirsts);
    free(eLasts);
    free(eThrones);
    free(timesDrawn);
    free(timesWon);
    free(timesDeadIfEWon);
    free(timesLivedIfEWon);
    free(timesThronedIfEWon);
    free(timesLastIfEWon);
}
