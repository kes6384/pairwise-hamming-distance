//Algorithm for finding the Hammond Distance profile of a given sequence
//Uses XOR and stores kmers in a hashtable
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "kc-c1.c" // k-mer counting
#include "khashl.h" // hash table

// Gets a sequence of length seqLen from the E. coli genome and encodes it as 2 bits for each character
// seqLen - length of sequence to get
// seq - unsigned int array to store bit-encoded sequence in
void getSequence(int seqLen , unsigned int* seq)
{
    FILE* genome = fopen("Escherichia_coli_0_1288_GCA_000303255.LargeContigs.fna" , "r");

    if(!genome)
        return;

    int charNum = 0; // Number of characters read so far
    while (charNum < seqLen)
    {
        int nxtChar = fgetc(genome);
        // check for end of file
        if(feof(genome))
            break;
        // check for header lines
        if(nxtChar == '>')
        {
            // move file pointer to next line
            while(nxtChar!='\n' && !(feof(genome)))
            {
                nxtChar = fgetc(genome);
            }
            if(feof(genome))
                break;
        }
        if(nxtChar != '\n')
        {
            int arrayIndex = charNum/16;
            if(nxtChar == 'A') // 00
                seq[arrayIndex] = seq[arrayIndex] | (0x00 << (charNum * 2));
            if(nxtChar == 'T') // 01
                seq[arrayIndex] = seq[arrayIndex] | (0x01 << (charNum * 2));
            if(nxtChar == 'C') // 10
                seq[arrayIndex] = seq[arrayIndex] | (0x02 << (charNum * 2));
            if(nxtChar == 'G') // 11
                seq[arrayIndex] = seq[arrayIndex] | (0x03 << (charNum * 2));
            charNum ++;
        }
    }

    fclose(genome);
}

// Helper function for hashing
// Calculates the key for a k-mer in the hash table
// Keys are based on the integer value of the bit-encoded k-mer
int getKey(unsigned int kmer , int kBits)
{
    // Largest possible key = kmer with all bits set
    unsigned int maxKey = (unsigned int)(pow(2 , kBits) - 1);

    return maxKey - kmer;
}

// Separate out and store all of the kmers in the sequence
// k-mers are subsequences of array of bits specified by seq
// k-mers start at bit index of seq specified by subseqStart param
// k-mers have length of kBits bits
// k-mers are stored in array kmers and are indexed by their starting position in the original sequence
// counts tracks the occurrences of each k-mer
void storeKmers(unsigned int *seq , int subseqStart , int kBits , unsigned int *kmers , unsigned int *counts)
{
    // Get subsequences, which are currently assumed to be <=16 characters (ie k <= 16)
    int subseqStartArrayPos = subseqStart/32;
    unsigned int subseq = (seq[subseqStartArrayPos] >> subseqStart);
    // Account for if the subsequence spans multiple array entries in seq
    if((kBits - (((subseqStartArrayPos + 1) * 32) - subseqStart)) > 0)
        subseq = subseq | (seq[subseqStartArrayPos + 1] << (((subseqStartArrayPos + 1) * 32) - subseqStart));

    // Cut subsequences down to correct length
    // Uses the fact that sum(2^n) from 0 to k is 2^(k)-1
    unsigned int cut = (unsigned int)(pow(2 , kBits) - 1);
    subseq = subseq & cut;

    counts[getKey(subseq , kBits)] += 1;

    // Store kmer
    kmers[subseqStart] = subseq;
}

// Find the Hamming distance of two k-mers of length k
// k-mers are subsequences stored with counts in kmers hash table
int hammingDist(uint32_t subseq1 , uint32_t subseq2 , int k)
{
    unsigned int result = subseq1 ^ subseq2;

    int dist = 0;
    for(int i=0; i<k; i++)
    {
        // Check character
        int nxtChar = result % 4;
        if(nxtChar != 0)
            dist ++;
        // Move to next character
        result = result >> 2;
    }

    return dist;
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// rate - sampling rate (theta1 * theta2)
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
void output(int *dists , int len)
{
    FILE* out = fopen("HD_XOR_output.txt" , "w");

    fprintf(out , "Hamming Distance : Number of Pairs");
    for(int i=0; i<=len; i++)
    {
        fprintf(out , "\n%d : %d" , i , (int)(dists[i]/(2)));
    }

    fclose(out);
}

// args: sequence length, k, file name
int main(int argc, char *argv[])
{
    int seqLen = atoi(argv[1]); // Sequence length
    int kVal = atoi(argv[2]); // k-mer length
    //int kBits = 2*kVal;
    // Hash table for storing all k-mer counts
    kc_c1_t *kmers;

    unsigned int *sequence = calloc(ceil(((float)seqLen*2)/8) , 1);
    char *file = argv[3];
    //char *file = "test.fna";
    kmers = count_file(file, kVal);
    getSequence(seqLen , sequence);

    // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
    int *dists = (int *)calloc(kVal+1 , sizeof(int));

    free(sequence);

    void **kmer1 = malloc(sizeof(void *));
    void **kmer2 = malloc(sizeof(void *));
    // Calculate Hamming distance for sampled k-mers
    for(int i = 0; i<kh_end(kmers); i +=1)
    {
        if(kh_exist(kmers, i))
        {
            //printf("!");
            uint32_t kmer1 = kh_key(kmers , i);
            //printf("%d\n" , kmer1);
            for(int j = 0; j<kh_end(kmers); j +=1)
            {
                if(kh_exist(kmers , j))
                {
                    uint32_t kmer2 = kh_key(kmers , j);
                    //printf("%d , %d\n" , kmer1 , kmer2);
                    if(kmer1 != kmer2)
                    {
                        int count1 = kh_val(kmers , kc_c1_get(kmers , kmer1));
                        int count2 = kh_val(kmers , kc_c1_get(kmers , kmer2));
                        dists[hammingDist(kmer1 , kmer2 , kVal)] += (count1 * count2);
                    }
                }
            }
        }
    }
    free(kmer1);
    free(kmer2);

    kc_c1_destroy(kmers);
    output(dists , kVal);
    free(dists);

    return 0;
}