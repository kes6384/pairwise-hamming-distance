//Naive algorithm for finding the Hammond Distance profile of a given sequence
//Takes FASTA files as input
//Outputs a text file representing the histogram of Hammond Distances between k-mers

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Gets a sequence of length seqLen from the E. coli genome
// seqLen - length of sequence to get
// seq - character array to store sequence in
void getSequence(int seqLen , char* seq , char* fileName)
{
    FILE* genome = fopen(fileName , "r");

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
        // append character to sequence array
        if(nxtChar != '\n')
        {
            seq[charNum] = nxtChar;
            charNum ++;
        }
    }

    fclose(genome);
}

// Find the Hamming distance of two strings of length len
int hammingDist(char *str1 , char *str2 , int len)
{
    int dist = 0;

    for(int i=0; i<len; i++)
    {
        if(str1[i] != str2[i])
        {
            dist ++;
        }
    }

    return dist;
}

// Output array of Hamming distance counts to a text file
// dists - array of Hamming distance counts
// len - length of array of Hamming distance counts (equivalent to k-mer size + 1, since HD ranges from 0 to k)
void output(int *dists , int len)
{
    FILE* out = fopen("naiveHammingDistanceOutput.txt" , "w");

    fprintf(out , "Hamming Distance : Number of Pairs");
    for(int i=0; i<=len; i++)
    {
        fprintf(out , "\n%d : %d" , i , dists[i]);
    }

    fclose(out);
}

// args: sequence length, k, file name
int main(int argc, char *argv[])
{
    int seqLen = atoi(argv[1]); // Sequence length
    int kVal = atoi(argv[2]); // k-mer length

    char *file = argv[3];
    char *sequence = calloc(seqLen , sizeof(char));
    getSequence(seqLen , sequence , file);

    // Tracks how many pairs had a Hamming distance of i, where i is an index of the array
    int *dists = (int *)calloc(kVal+1 , sizeof(int));

    // Calculate Hamming distance for all k-mer pairs in the sequence
    for(int i = 0; i < (seqLen-kVal)+1; i +=1)
    {
        for(int j = i+1; j < (seqLen-kVal)+1; j +=1)
        {
            if(i!=j)
            {
                dists[hammingDist(&sequence[i] , &sequence[j] , kVal)] ++;
            }
        }
    }

    // Ignore duplicates
    dists[0] = 0;

    free(sequence);
    output(dists , kVal);
    free(dists);

    return 0;
}