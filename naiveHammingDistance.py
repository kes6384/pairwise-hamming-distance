#Naive algorithm for finding the Hammond Distance profile of a given sequence
#Kathryn Silverberg

#Takes FASTA files as input
#Outputs a text file representing the histogram of Hammond Distances between k-mers

#Gets a sequence of length seqLen from the E. coli genome
def getSequence(seqLen):
    file = open("Escherichia_coli_0_1288_GCA_000303255.LargeContigs.fna" , "r")

    seq = ""
    while len(seq) < seqLen:
        nextLine = file.readline().replace("\n" , "")
        if ">" not in nextLine:
            seq += nextLine

    file.close()
    return seq

#Find the Hamming distance of two strings
def hammingDist(str1 , str2):
    dist = 0
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            dist += 1
    return dist

#Output the Hamming distance counts to a text file
def output(dists):
    out = open("naiveHammingDistanceOutputPYTHON.txt" , "w")
    out.write("Hamming Distance : Number of Pairs")
    for i in range(0,len(dists)):
        out.write(f"\n{i} : {dists[i]}")
    out.close()

def main():
    seqLen = 10000 #Sequence length
    kVal = 16 #k-mer length

    sequence = getSequence(seqLen)
    #Tracks how many pairs had a Hamming distance of i, where i is the list index
    dists = [0 for i in range(0,kVal+1)]

    #Calculate HD for all k-mer pairs in the sequence
    for i in range(0,seqLen-kVal+1):
        for j in range(i+1,seqLen-kVal+1):
            if(i!=j):
                dists[hammingDist(sequence[i:i+kVal],sequence[j:j+kVal])] += 1

    output(dists)

main()
