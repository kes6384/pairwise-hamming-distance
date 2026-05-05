# Hamming Distance Profile

There are two main options offered for calculating the pairwise Hamming distance of a DNA sequence.
Sequences are assumed to be in FASTA file format.

## Sampling Method

The sampling option allows for users to compute either the full exact pairwise Hamming distance profile or an estimate. The source file for the sampling method is bitEncoding_sampling.c, and once compiled can be called via the command line. Command line execution should include arguments for k, theta1, theta2, and the FASTA file. Theta1 and theta2 represent the sampling rates, and should be at most 1.0. If both are set to 1.0, then the result of running the command will be the exact Hamming distance profile. Otherwise, theta1 and theta2 represent the rate at which each k-mer in a pair are sampled for the computation of Hamming distance. Results will be output to a text file with the name "HDSampling_output.txt" and will be in the format of two columns (Hamming Distance : Number of Pairs).

## Hamming Distance of 1

For users who are interested specifically in the number of k-mer pairs with a Hamming distance of 1, a more efficient option is offered. This method will result in counting only the number of pairs with a Hamming distance of 1, thus allowing for much faster results. The source file for this method is HD1.c, and once compiled can be called via the command line. Command line execution should include arguments for k and the FASTA file. Results will be output to a text file with the name "HD1_Output.txt" and will be in the format of two columns (Hamming Distance : Number of Pairs).

## Building Files

bitEncoding_sampling.c and HD1.c should both be compiled using the following command: gcc -O2 -Wno-unused-function -o [EXECUTABLE NAME] [SOURCE FILE NAME].c -lm -lz -lpthread. For the sampling method, the additional files murmurhash.c, murmurhash.h, kc-c1.c, and khash1.h files are required. For the Hamming distance of 1 method, only the kc-c1.c and khash1.h files are required.
