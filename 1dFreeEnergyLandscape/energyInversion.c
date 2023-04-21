#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#define R1 8.31446261815324 // J.K^-1.mol^-1 - GROMACS
#define R2 1.98720425864083 // cal.K^-1.mol^-1 - LAMMPS

typedef struct probabilityDistribution
{
	int nHbonds;
	float probability, energy;
} PROBABILITY_DISTRIBUTION;

PROBABILITY_DISTRIBUTION *initializeDistData (PROBABILITY_DISTRIBUTION *distData, int nLines, int nHbondsMax, int nHbondInterval)
{
	for (int i = 0; i < nLines; ++i)
	{
		distData[i].probability = 0;
		distData[i].energy = 0;

		if (i == 0) {
			distData[i].nHbonds = 0; }
		else {
			distData[i].nHbonds = distData[i - 1].nHbonds + nHbondInterval; }
	}

	return distData;
}

PROBABILITY_DISTRIBUTION *readData (FILE *input, PROBABILITY_DISTRIBUTION *distData, int nHbondsMax, int nHbondInterval, int nLines)
{
	rewind (input);

	char lineString[2000];
	int tempX, index;
	float tempY;

	while (fgets (lineString, 2000, input) != NULL)
	{
		sscanf (lineString, "%d %f\n", &tempX, &tempY);

		if (tempX == 0)
		{
			distData[0].probability += tempY;
		}
		else
		{
			index = tempX / nHbondInterval;
			distData[index].probability += tempY;
		}
	}

	return distData;
}

PROBABILITY_DISTRIBUTION *computeEnergy (PROBABILITY_DISTRIBUTION *distData, int nLines, float temperature)
{
	for (int i = 0; i < nLines; ++i)
	{
		if (distData[i].probability == 0)
		{
			distData[i].energy = 0;
		}
		else
		{
			distData[i].energy = R1 * temperature * log (distData[i].probability);
		}
	}

	return distData;
}

void printDistData (PROBABILITY_DISTRIBUTION *distData, int nLines, FILE *output)
{
	for (int i = 0; i < nLines; ++i)
	{
		fprintf(output, "%d %f %f\n", distData[i].nHbonds, distData[i].probability, distData[i].energy);
	}
}

int main(int argc, char const *argv[])
{
	/*
	REQUIRED ARGS:
	~~~~~~~~~~~~~

	argv[0] = program
	argv[1] = input file name (*.dat)
	argv[2] = output file name
	argv[3] = Max. H bonds
	argv[4] = hBonds interval
	argv[5] = temperature (in Kelvins)

	*/
	FILE *input, *output;
	input = fopen (argv[1], "r");
	output = fopen (argv[2], "w");

	int nHbondsMax = atoi (argv[3]), nHbondInterval = atoi (argv[4]);
	float temperature = atof (argv[5]);

	int nLines = (nHbondsMax / nHbondInterval) + 1;

	PROBABILITY_DISTRIBUTION *distData;
	distData = (PROBABILITY_DISTRIBUTION *) malloc (nLines * sizeof (PROBABILITY_DISTRIBUTION));

	distData = initializeDistData (distData, nLines, nHbondsMax, nHbondInterval);
	distData = readData (input, distData, nHbondsMax, nHbondInterval, nLines);

	distData = computeEnergy (distData, nLines, temperature);
	printDistData (distData, nLines, output);

	free (distData);
	fclose (input);
	fclose (output);
	return 0;
}