# LGM-attack

C code implementing a key recovery attack on the scheme by Li, Galbraith and Ma for fully homomorphic encryption. The [LGM scheme](https://eprint.iacr.org/2016/1146.pdf) uses a list of vectors drawn from a Gaussian distribution as the secret key. The attack tries to recover coefficients from these vectors by asking for many decryptions of the same ciphertexts and using statistical methods to remove the uncertainty introduced by the random elements used in the decryption function.

## Using the Code

### generateEvectors.c

The program generateEvectors.c is used to sample a number of vectors used as the secret key. This program requires the [Discrete Gaussian Sampler]( https://bitbucket.org/malb/dgs/src/master/) library to work. Running the code takes 3 parameters:

* Length of each vector
* Number of vectors produced
* Standard deviation for the discrete Gaussian sampler

The output is a file containg the secret key vectors that can be read by coefficientSearch.c

### coefficientSearch.c

The is the main program that simulates LGM decryption and runs the attack trying to recover coefficients from the secret vectors. The program needs a file with the secret key vectors to run, but otherwise only depends on standard C libraries.  The program takes six parameters to run:

* lambda_max: the maximum value a random lambda_i can take in the decryption function (typically 2 or 3)

* Sample size: the number of decryption queries used for each ciphertext to generate a count of the number of 1-decryptions. Larger values give more accurate results, depending on the length of the secret vectors and the standard deviation used for sampling them.

* negative: flag that indicates the interval for sampling random lambda-values.  0 means the lambda_i are sampled from [0...lambda_max-1] and 1 means the lambda_i are sampled from [-lambda_max...lambda_max].

* start: which vector in the file to start searching coefficients for

* stop: when to stop processing vectors from the file. Start and stop allow for easy parallellization.

* file: filename of file with secret key vectors.

The program outputs a file with the estimated coefficients, and notes the difference with the correct value if the estimate was not correct.

