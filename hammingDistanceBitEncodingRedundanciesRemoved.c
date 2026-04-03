//Naive algorithm for finding the Hammond Distance profile of a given sequence
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

// Separate out and store all of the kmers in the sequence
// k-mers are subsequences of given array of bits specified by seq
// k-mers start at bits specified by subseqStart param
// k-mers have length of kBits bits
// k-mers are stored in the array kmers and are indexed by their starting position in the original sequence
void storeKmers(unsigned int *seq , int subseqStart , int kBits , unsigned int *kmers)
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

    // Store kmer
    kmers[subseqStart] = subseq;
}

// Find the Hamming distance of two k-mers of length kBits bits
// k-mers are subsequences of given array of bits stored in kmers array
// k-mers start at bits specified by subseq1Start and subseq2Start params
int hammingDist(unsigned int *kmers , int subseq1Start , int subseq2Start , int kBits)
{
    // Get subsequences out of kmers array
    unsigned int subseq1 = kmers[subseq1Start];
    unsigned int subseq2 = kmers[subseq2Start];

    // Mask bits to separate first and second bits of each character for each subsequence
    int bitMaskLen = ceil(((float)kBits)/8);
    unsigned int *mask = calloc(bitMaskLen , 1);
    memset(mask , 170 , bitMaskLen);
    unsigned int subseq1_firstBit = subseq1 & *mask;
    unsigned int subseq2_firstBit = subseq2 & *mask;
    unsigned int subseq1_secondBit = subseq1 & (*mask >> 1);
    unsigned int subseq2_secondBit = subseq2 & (*mask >> 1);

    // XOR each masked subsequence to compare first and second bits
    // XOR result is 0 if the bits match
    unsigned int xor_firstBits = subseq1_firstBit ^ subseq2_firstBit;
    unsigned int xor_secondBits = (subseq1_secondBit ^ subseq2_secondBit) << 1;

    // NOR = 1 if both bits for the character matched = the character matched
    // ~NOR = 1 if the characters did not match
    // popcount(~nor) = number of characters that did not match
    // ~nor = or
    unsigned int or = xor_firstBits | xor_secondBits;
    int dist = __builtin_popcount(or);

    return dist;
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
void output(int *dists , int len)
{
    FILE* out = fopen("bitEncodingHammingDistanceOutput.txt" , "w");

    fprintf(out , "Hamming Distance : Number of Pairs");
    for(int i=0; i<=len; i++)
    {
        fprintf(out , "\n%d : %d" , i , dists[i]);
    }

    fclose(out);
}

void main()
{
    int seqLen = 10000; // Sequence length
    int kVal = 16; // k-mer length

    unsigned int *sequence = calloc(ceil(((float)seqLen*2)/8) , 1);
    getSequence(seqLen , sequence);

    // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
    int *dists = (int *)calloc(kVal+1 , sizeof(int));

    // Separate out all k-mers
    int kBits = 2*kVal;
    unsigned int *kmers = calloc(seqLen-kVal , kBits);
    for(int i = 0; i < (seqLen-kVal)*2+1; i +=2)
    {
        storeKmers(sequence , i , kBits , kmers);
    }

    // Calculate Hamming distance for all k-mer pairs in the sequence
    for(int i = 0; i < (seqLen-kVal)*2+1; i +=2)
    {
        for(int j = i+2; j < (seqLen-kVal)*2+1; j +=2)
        {
            if(i!=j)
            {
                dists[hammingDist(kmers , i , j , kBits)] ++;
            }
        }
    }

    free(sequence);
    free(kmers);
    output(dists , kVal);
    free(dists);
}