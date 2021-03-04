# LGMattack

C-code implementing an attack on the LGM scheme for fully homomorphic encryption.  The LGM scheme uses a list of vectors drawn from a Gaussian distribution as the secret key.  The attack tries to recover coefficients from these vectors by asking for many decryptions of the same ciphertexts and using statistical methods to remove the uncertainty introduced by the random elements used in the decryption function.

## Using the Code

### generateEvectors.c

The program generateEvectors.c is used to sample a number of vectors used as the secret key.  This program depends itself on the Discrete Gaussian Sampler at https://bitbucket.org/malb/dgs/src/master/ to work.  Running the code takes 3 parameters:

* Length of each vector
* Number of vectors produced
* Standard deviation for the Gaussian sampler

The output is a file containg the vectors that can be read by coefficientSearch.c

### coefficientSearch.c

The is the main program that simulates LGM decryption and runs the attack trying to recover coefficients from the secret vectors.  The program needs a file with the secret vectors to run, but otherwise only depends on standard C libraries.  The program takes six parameters to run:

* lambda-max: the maximum value a random lambda-i can take in the decryption function (typically 2 or 3)

* Sample size: the number of decryption queries used for each ciphertext to generate a count of the number of 1-decryptions.  Larger values give more accurate results, depending on the length of the secret vectors and the standard deviation used for sampling them.

* negative: flag that indicates the interval for sampling random lambda-values.  0 means lambda-i sampled from [0...lambda-max-1] and 1 means lambda-i sampled from [-lambda-max...lambda-max].

* start: which vector in the file to start searching coefficients for

* stop: when to stop processing vectors from the file.  Start and stop allow for easy parallellization.

* file: filename of file with secret vectors.

The program outputs a file with the estimated coefficients, and notes the difference with the correct value if the estimate was not correct.

