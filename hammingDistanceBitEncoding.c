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
            if(nxtChar == 'A') // 00
                seq[charNum/16] = seq[charNum/16] | (0x00 << (charNum * 2));
            if(nxtChar == 'T') // 01
                seq[charNum/16] = seq[charNum/16] | (0x01 << (charNum * 2));
            if(nxtChar == 'C') // 10
                seq[charNum/16] = seq[charNum/16] | (0x02 << (charNum * 2));
            if(nxtChar == 'G') // 11
                seq[charNum/16] = seq[charNum/16] | (0x03 << (charNum * 2));
            charNum ++;
        }
    }

    fclose(genome);
}

// Find the Hamming distance of two k-mers of length k
// k-mers are subsequences of given array of bits specified by seq
// Sequence has length len
// k-mers start at bits specified by subseq1Start and subseq2Start params
int hammingDist(unsigned int *seq , int len , int subseq1Start , int subseq2Start , int k)
{
    // Get subsequences, which are currently assumed to be <=16 characters (ie k <= 16)
    unsigned int subseq1 = (seq[subseq1Start/32] >> subseq1Start);
    unsigned int subseq2 = (seq[subseq2Start/32] >> subseq2Start);
    // Account for if the subsequence spans multiple array entries in seq
    if((2*k - (((subseq1Start/32 + 1) * 32) - subseq1Start)) > 0)
        subseq1 = subseq1 | (seq[subseq1Start/32 + 1] << (((subseq1Start/32 + 1) * 32) - subseq1Start));
    if((2*k - (((subseq2Start/32 + 1) * 32) - subseq2Start)) > 0)
    {
        subseq2 = subseq2 | (seq[subseq2Start/32 + 1] << (((subseq2Start/32 + 1) * 32) - subseq2Start));
    }

    // Cut subsequences down to correct length
    // Uses the fact that sum(2^n) from 0 to k is 2^(k)-1
    subseq1 = subseq1 & (unsigned int)(pow(2 , k*2) - 1);
    subseq2 = subseq2 & (unsigned int)(pow(2 , k*2) - 1);

    // Mask bits to separate first and second bits of each character for each subsequence
    unsigned int *maskFirstBit = calloc(ceil(((float)k*2)/8) , 1);
    unsigned int *maskSecondBit = calloc(ceil(((float)k*2)/8) , 1);
    memset(maskFirstBit , 170 , ceil(((float)k*2)/8));
    memset(maskSecondBit , 85 , ceil(((float)k*2)/8));
    unsigned int subseq1_firstBit = subseq1 & *maskFirstBit;
    unsigned int subseq2_firstBit = subseq2 & *maskFirstBit;
    unsigned int subseq1_secondBit = subseq1 & *maskSecondBit;
    unsigned int subseq2_secondBit = subseq2 & *maskSecondBit;

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

    // Calculate Hamming distance for all k-mer pairs in the sequence
    for(int i = 0; i < (seqLen-kVal)*2+1; i +=2)
    {
        for(int j = i+2; j < (seqLen-kVal)*2+1; j +=2)
        {
            if(i!=j)
            {
                dists[hammingDist(sequence , seqLen , i , j , kVal)] ++;
            }
        }
    }

    free(sequence);
    output(dists , kVal);
    free(dists);
}