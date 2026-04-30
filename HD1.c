//Algorithm for finding the number of k-mer pairs in a given sequence with Hamming Distance of 1
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "kc-c1.c"
#include "khashl.h" // hash table
KHASHL_MAP_INIT(, kc_c1_t, kc_c1, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

// Gets a sequence of length seqLen from the E. coli genome and encodes it as 2 bits for each character
// seqLen - length of sequence to get
// seq - unsigned int array to store bit-encoded sequence in
// k - kmer length
void getSequence(int seqLen , unsigned int* seq , int k , kc_c1_t *kmers)
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
    // Store all k-mers in hash table and find kmer counts
    kmers = count_file(genome, k);
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

// Check if there are any pairs of kmers that include paramter kmer with HD of 1
int hammingDist(unsigned int kmer , int kBits , kc_c1_t* kmers)
{
    int count = 0;
    // Try swapping out one character in the k-mer to another one in the alphabet
    // Check if result is in the kmer hash table
    unsigned int partner = 0;
    for(int i = 0; i < (kBits/2); i+=2)
    {
        // Build partner by replacing one location of k-mer with other possible characters
        partner = 0;
        partner = partner | ((kmer >> (i+2)) << (i+2));
        partner = partner | 00 << (i);
        if(i != 0)
            partner = partner | (kmer % (4*i));
        if(partner != kmer)
        {
            if(prefix_get(kmers, partner) == kh_end(kmers))
                count++;
        }
        partner = 0;
        partner = partner | ((kmer >> (i+2)) << (i+2));
        partner = partner | 01 << (i);
        if(i != 0)
            partner = partner | (kmer % (4*i));
        if(partner != kmer)
        {
            if(prefix_get(kmers, partner) == kh_end(kmers))
                count++;
        }
        partner = 0;
        partner = partner | ((kmer >> (i+2)) << (i+2));
        partner = partner | 10 << (i);
        if(i != 0)
            partner = partner | (kmer % (4*i));
        if(partner != kmer)
        {
            if(prefix_get(kmers, partner) == kh_end(kmers))
                count++;
        }
        partner = 0;
        partner = partner | ((kmer >> (i+2)) << (i+2));
        partner = partner | 11 << (i);
        if(i != 0)
            partner = partner | (kmer % (4*i));
        if(partner != kmer)
        {
            if(prefix_get(kmers, partner) == kh_end(kmers))
                count++;
        }
    }

    return count;
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
void output(int count , int len)
{
    FILE* out = fopen("bitEncodingHD1_Output.txt" , "w");

    fprintf(out , "Hamming Distance : Number of Pairs");
    fprintf(out , "\n1 : %d" , count);

    fclose(out);
}

void main()
{
    int seqLen = 20; // Sequence length
    int kVal = 2; // k-mer length
    int kBits = 2*kVal;
    // Hash table for storing k-mers
    kc_c1_t *kmers;

    // Get sequence and store kmers
    unsigned int *sequence = calloc(ceil(((float)seqLen*2)/8) , 1);
    getSequence(seqLen , sequence , kVal , kmers);

    free(sequence);

    // Tracks how many pairs had a Hamming distance of 1
    int count = 0;
    // Check each k-mer for HD = 1 partners
    khint_t pos;
    for(int pos=0; pos<kh_end(kmers); pos++)
    {
        uint64_t kmer = kh_key(kmers, pos);
        count += hammingDist(kmer , kBits , kmers);
    }

    kc_c1_destroy(kmers);
    output(count/2 , kVal);
}