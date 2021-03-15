/**
 * \file coefficientSearch.c
 *
 * \author H{\aa}vard Raddum, Prastudy Fauzi
 *
 * The is the main program that simulates LGM decryption and runs the key recovery attack.
 * The attack tries to recover individual coefficients e_{i,j} from the secret key vectors e_i.
 * The program needs a file containing secret key vectors to run, but otherwise only depends on standard C libraries
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

/* In case random() and srandom() not available (e.g, in Windows) */
#ifndef srandom
#if defined(__MINGW32__) || defined(__MINGW64__) || defined(WIN32)
#warning "Using rand and srand instead of random and srandom"
#define random() rand()
#define srandom(unsigned) srand(unsigned)
#endif
#endif

int sampleSize; // the number of decryption queries T required for the attack to succeed
int lambda_max; // determines the maximum value b of the random vectors lambda_i used to generate a one-time secret key 
int negative; // 0 if lambda_i's are from {0,1,...,b-1} and 1 if lambda_i's are from {-b,...,b}
int q = 94980001; // main parameter q, the prime modulus
int qBy4; // essentially 2^{j-2}, where j is the length of q in bits
double Atilde; // key coefficient required to estimate the e_{i, j}'s


int *randomVector(int ne){
	/* Random non-zero vector of size ne and values in [0...lambda_max-1] or [-lambda_max ... lambda_max] */
	int * v = (int *) calloc(ne, sizeof(int));
	int width = lambda_max; // default: valid values [0...lambda_max-1]
	int min = 0;
	
	if(negative){ // valid values are [-lambda_max...lambda_max]
		width = 2 * lambda_max + 1; 
		min = - lambda_max;
	}
	
	for(int i=0; i < ne; i++){
		v[i] = (random() % width) + min;
	}
	
	// Check if we get a zero vector
	int zero = 1;
	for(int i = 0; i < ne; i++){
		if(v[i]){
			zero = 0;
			break;
		}
	}
	
	// If we get a zero vector, set a random position to be 1
	if(zero){
		int i = random() % ne;
		v[i] = 1; 
	}
	
	return v;
}

double norm(int *e, int m){
	// Returns the Euclidean norm of e. */
	int i, sum = 0;

	for(i = 0; i < m; i++){
		sum += (e[i] * e[i]);
	}
	
	return sqrt(sum);
}

int innerProduct(int *e, int n, int *b){
	/* Computes and returns <e,b>, where b is a vector of length n. */
	int i, sum = 0;
	
	for(i = 0; i < n; i++){
		sum += (e[i] * b[i]);
	}

	return sum;
}

int decryptRowMatrix(int *e, int ne, int attacked_coefficient, int ai){
	/*	Selects lambda-vector at random, and returns the decryption result of the ciphertext with ai's in row of attacked_coefficient, and 1's in row t+1.
			Returns 0 if lambda_{attacked_coefficient} = 0 or |<e, lambda> + ai| < qBy4,
			and returns 1 if lambda_{attacked_coefficient} = 1 and |<e, lambda> + ai| >= qBy4
	*/
	int *lambda;
	int u;
	
	lambda = randomVector(ne);
	u = lambda[attacked_coefficient] * ai + innerProduct(e, ne, lambda);
		
	free(lambda);
	if(u < qBy4 && u > - qBy4){ // |u| < cutoff
		return 0;
	}

	return 1;
}

int decryptDiagonalMatrix(int *e, int ne, int ai){
	/*	Selects lambda-vector at random, and returns the decryption result of the ciphertext matrix with ai's in position (i,i) of the matrix, and 1's in row t + 1.
			Will always give u = ai + <e, lambda>, regardless of the lambda or chosen column.
			Returns 0 if |u| < qBy4, and 1 otherwise.
	*/
	int *lambda;
	double u;
	
	lambda = randomVector(ne);
	int i = random() % ne;
	
	while(!lambda[i]) { // ensure we get a random non-zero element lambda[i]
		i = (i + 1) % ne;
	}
	
	u = lambda[i] * ai + innerProduct(e, ne, lambda);
	free(lambda);
	
	if(u < qBy4 && u > - qBy4) // |u| < cutoff
		return 0;

	return 1;    
}

int numberOf1Decryptions(int *e, int ne, int attacked_coefficient, int ai, const char *matrix_type){
	/* Number of times decryption oracle returns 1 out of sample_size queries. */
	int n1 = 0, i;

	if(strcmp(matrix_type, "row_ai") == 0){
		for(i = 0; i < sampleSize; i++){
			if(decryptRowMatrix(e, ne, attacked_coefficient, ai))
				n1++;
		}
	}
	
	if(strcmp(matrix_type, "diagonal_ai") == 0){
		for(i = 0; i < sampleSize; i++){
			if(decryptDiagonalMatrix(e, ne, ai))
				n1++;
		}
	}

	return n1;
}


double search_a(int *e, int ne, int attacked_coefficient, double scalingFactor, const char *matrixType){
	/* Using binary search to estimate \tilde{a} or \tilde{A} */
	int min_a, max_a, middle;
	int n1, nmin, nmax, Nhalf;
	int jump;
	double bias;
	int middle_a;
 
	Nhalf = sampleSize/2;

	if(strcmp(matrixType, "diagonal_ai") == 0){ // Searching for Atilde, use qBy4/lambda_max as starting point for binary search
		if(negative)
			middle_a = qBy4/lambda_max;
		else
			middle_a = qBy4/(lambda_max-1);
	}
	else // searching for atilde, use Atilde as starting point for binary search
		middle_a = Atilde;

	n1 = scalingFactor * numberOf1Decryptions(e, ne, attacked_coefficient, middle_a, matrixType);
	jump = 1;
	if(n1 < Nhalf){ //middle_a gives n1 larger than N/2, jumping upwards with increasingly large jumps until nmax>Nhalf
		min_a = middle_a;
		nmin = n1;
		max_a = min_a + jump;
		nmax = scalingFactor * numberOf1Decryptions(e, ne, attacked_coefficient, max_a, matrixType);
		while(nmax<Nhalf){
			min_a = max_a;
			nmin = nmax;
			jump *= 2; // double the size of the next jump
			max_a = max_a + jump;
			nmax = scalingFactor * numberOf1Decryptions(e, ne, attacked_coefficient, max_a, matrixType);
		} // here we know that nmin < N/2 <= nmax, after a logarithmic number of iterations in the while loop
	}
	else{ // middle_a gives n1 smaller than N/2, jumping downwards until nmin < Nhalf
		max_a = middle_a;
		nmax = n1;
		min_a = max_a - jump;
		nmin = scalingFactor * numberOf1Decryptions(e, ne, attacked_coefficient, min_a, matrixType);
		while(nmin>Nhalf){
			max_a = min_a;
			nmax = nmin;
			jump *= 2; // double the size of the next jump
			min_a = min_a - jump;
			nmin = scalingFactor * numberOf1Decryptions(e, ne, attacked_coefficient, min_a, matrixType);
		} // here we know that nmin <= N/2 < nmax, after a logarithmic number of iterations in the while loop
	}
	
	// Starting actual binary search
	while(min_a < max_a - 1) {
		middle = (min_a + max_a)/2;
		n1 = scalingFactor * numberOf1Decryptions(e, ne, attacked_coefficient, middle, matrixType);
		if(n1 > Nhalf){
			max_a = middle;
			nmax = n1;
		}
		else{
			min_a = middle;
			nmin = n1;
		}
	} // here nmin < N/2 <= nmax, and min_a = max_a - 1
	
	// Interpolation
	bias = ((double)(nmax - Nhalf)) / ((double)(nmax - nmin));
	
	return bias * (double) min_a + (1.0 - bias) * (double) max_a;
}

void checkE(char *destFileName, int **Etab, int *all_e, int attacked_vector, int m){
	/* 	Compares estimated vector all_e with the correct secret key Etab[.][attacked_vector].
			The function outputs these vectors, the number of incorrect positions and the Euclidean norm of their distance. 
	*/
	int wrong = 0;
	int *diff_vector = (int *) malloc(m * sizeof(int));
	FILE *fp = fopen(destFileName, "a");
	fprintf(fp, "output s_%d = [", attacked_vector);
	for(int i = 0; i < m; i++){
		if(all_e[i] == (double) Etab[i][attacked_vector])
			fprintf(fp, "%d ", all_e[i]);
		else {
			fprintf(fp, "%d* ", all_e[i]);
			wrong++;
		}
		diff_vector[i] = all_e[i] - Etab[i][attacked_vector];
	}
	fprintf(fp, "]\n");
	
	fprintf(fp, "actual s_%d = [", attacked_vector);
	for(int i = 0; i < m; i++){
		fprintf(fp,"%d ",Etab[i][attacked_vector]);
	}
	fprintf(fp, "]\n");
	fprintf(fp, "Wrong: %d positions\n", wrong);
	fprintf(fp, "Distance from real s_%d: %1.5f \n", attacked_vector, norm(diff_vector, m));
	fclose(fp);
}
	
int main(int argc, char *argv[]){
	/**
	*	Program parameters:
	*	- lambda_max (e.g., 2)
	*	- sample size (e.g. 100000)
	*	- 0 (if lambda in [0...lambda_max - 1]) or 1 if lambda in [0...lambda_max])
	*	- which of the secret key vectors e_i to recover, i in [0...t-1] lambda_max
	* - start and stop of which coefficients to find (like 24 27 for coefficients 24, 25 and 26),
	* - the secret key file (e.g., secretKey_100_75_25.txt)
	* 
	* The result will be written to a file with the format "estCoeff%dmill_%d_%d_%d_s%d_b%d.txt"
	* For example, if recovering e_3 with the example parameters above, the file name will be "estCoeff0mill_100_75_25_s3_b2.txt"
	*/
	FILE *fp;
	int i, j, qj = ceil(log2(q));
	int **Etab, attacked_vector, attacked_index, t, m, sum, sigma;
	int start, stop; // which coefficients of secret key vectors to estimate

	qBy4 = 1 << (qj - 2); // 2^{j-2}
	srandom(time(NULL));
	
	if(argc!=8){
		printf("Program must be called with:\n - lambda_max as first parameter (like 2, 3 or 4), \n - sample size as second parameter (like 100000, 0 for default sample size),\n - 0 (for lambda in [0...b-1]) or 1 (for lambda in [-b...b]) as third parameter,\n - which of the secret key vectors s_i to recover, i in [0...t-1],\n - start and stop of which coefficients to find (like 24 27 for coefficients 24, 25 and 26),\n the secret key file (e.g., secretKey_100_75_25.txt).\n");
		exit(0);
	}

	negative = 0;
	
	lambda_max=  atoi(argv[1]);
	sampleSize = atoi(argv[2]);
	negative = atoi(argv[3]);
	attacked_vector = atoi(argv[4]);
	start = atoi(argv[5]);
	stop = atoi(argv[6]);

	sscanf(argv[7], "secretKey_%*d_%*d_%d", &sigma);
	fp = fopen(argv[7], "r");
	fscanf(fp, "m=%d t=%d\n", &m, &t);
	
	for(i = 0; i < start; i++){ // ignore secret key coefficients that won't be attacked 
		fscanf(fp, "[");
		for(j = 0; j < t; j++)
			fscanf(fp, "%*d ");
		fscanf(fp, "]\n");
	}
	
	Etab = (int **) malloc((stop-start) * sizeof(int *));
	for(i = 0; i < stop-start; i++){
		Etab[i] = (int *) malloc(t * sizeof(int)); // i-th element of each s_j
		fscanf(fp, "[");
		for(j = 0; j < t; j++)
			fscanf(fp, "%d ", Etab[i] + j);
		fscanf(fp, "]\n");
	}
	fclose(fp);
	
	double recovered_e, estimated_e, scaleDiag, scaleRow, atilde;

	if(negative){
		scaleDiag = (double) lambda_max;
		scaleRow = ((double) (2 * lambda_max + 1)) / 2.0;
	}
	else{
		scaleDiag = (double) (lambda_max - 1);
		scaleRow = (double) lambda_max;
	}
	
	int* all_e = (int *) malloc(m * sizeof(int));
	char filename[80];
	sprintf(filename, "estCoeff%dmill_%d_%d_%d_s%d_b%d.txt", sampleSize/1000000, t, m, sigma, attacked_vector, lambda_max);
	
	
	for(attacked_index = 0; attacked_index < stop - start; attacked_index++){ // implement the key recovery algorithm described in Section 4 of our paper
		Atilde = search_a(Etab[attacked_index], t, attacked_vector, scaleDiag, "diagonal_ai"); // default: Atilde = 2^{j-2} - sum(e_k)/2 + small error

		atilde = search_a(Etab[attacked_index], t, attacked_vector, scaleRow, "row_ai"); // default: atilde = 2^{j-2} - sum(e_k)/2 - e_i/2 + small error
		
		if(!negative)
			estimated_e = 2 * (Atilde - atilde);
		else
			estimated_e = Atilde - atilde;
		
		fp = fopen(filename, "a");
		fprintf(fp, "estimated e_%d,%d = %1.3f", attacked_vector, attacked_index + start, estimated_e);
		recovered_e = round(estimated_e); // gives the correct coefficient assuming the errors in Atilde and atilde are small
		all_e[attacked_index] = recovered_e;
		if(recovered_e == (double) Etab[attacked_index][attacked_vector])
			fprintf(fp, "\n");
		else
			fprintf(fp, ", off by: %1.1f\n", recovered_e - (double) Etab[attacked_index][attacked_vector]);
		fclose(fp);
	}
	if(stop - start == m)
		checkE(filename, Etab, all_e, attacked_vector, m); // print the whole recovered e vector and the actual secret key vector as comparison 
}
	
