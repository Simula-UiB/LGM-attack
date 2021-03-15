/**
 * \file generateEvectors.c
 *
 * \author H{\aa}vard Raddum, Prastudy Fauzi
 *
 * Code for sampling LGM secret keys.
 * It requires the DGS library https://github.com/malb/dgs to generate discrete Gaussians over the integers 
 */

#include <dgs/dgs_gauss.h>

int *randomE(int ne, double sigma, int B){
	/* Samples a vector of integers from a Gaussian distribution with standard deviation sigma and cut-off bound B. */
	dgs_disc_gauss_dp_t *DGS;
	int *e, i;

	e = (int *) malloc(ne * sizeof(int));
	DGS = dgs_disc_gauss_dp_init(sigma, 0.0, (size_t) B, DGS_DISC_GAUSS_DEFAULT);

	for(i = 0; i < ne; i++)
		e[i] = (int) DGS->call(DGS);

	dgs_disc_gauss_dp_clear(DGS);

	return e;
}

int main(int argc, char *argv[]){
	/**
	* Program parameters:
	* - number of secret key vectors t (e.g., 190)
	* - length m of secret ley vectors (e.g., 525)
	* - the standard deviation (sigma) for the discrete Gaussian sampler (e.g., 25)
	* 
	* The result will be written to a file with the format "secretKey_%d_%d_%d.txt"
	* For example, if t = 190, m = 525, sigma = 25 the file name will be "secretKey_190_525_25.txt"
	*/	
	FILE *fp;
	char filename[80];
	int **Etab;
	int t, m, sigma, i, j;
	int B = 150;

	if(argc<4){
		printf("Usage: specify number of vectors produced (t) as first parameter, length of each vector (m) as second parameter, and standard deviation (sigma) for the Gaussian sampler as the last parameter.\n"); exit(0);
	}

	t = atoi(argv[1]);
	m = atoi(argv[2]);
	sigma = (double) atoi(argv[3]);

	sprintf(filename, "secretKey_%d_%d_%d.txt", t, m, sigma);
	fp = fopen(filename, "w");
	fprintf(fp, "m=%d t=%d\n", m, t);

	Etab=(int **) malloc(t * sizeof(int *));
	for(i = 0; i < t; i++)
		Etab[i] = randomE(m, (double) sigma, B);

	/*
	Vectors are stored as m vectors of length t in the file for easy parallelization.
	The first vector contains [e_{1,0} e_{2,0} ... e_{t,0}], the second vector contains [e_{1,1} e_{2,1} ... e_{t,1}], etc.
	*/
	for(j = 0; j < m; j++){
		fprintf(fp, "[");
		for(i = 0; i < t; i++)
			fprintf(fp, "%d ", Etab[i][j]);
		fprintf(fp, "]\n");
	}
	fclose(fp);
}
