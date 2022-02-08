import math
import random
from Bio import Entrez, SeqIO

def getCounts(motifs):
    """Return nt counts given a list of motif strings."""
    
    counts = []
    for i in range(len(motifs[0])):
        counts.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
    
    for motif in motifs:
        motif = motif.upper()
        for i in range(len(motif)):
            counts[i][motif[i]] += 1
    
    return counts
    
def getProfile(motifs):
    """Get a profile from a set of motifs."""
    
    counts = getCounts(motifs)
    total = len(motifs)
    for i in range(len(counts)):
        for nt in 'ACGT':
            counts[i][nt] /= total
    return counts  # really now a profile
    
def getProfileLaplace(motifs):
    """Get a profile from a set of motifs, adding pseudocounts 
       (Laplace's rule of succession) to prevent zero-probability events.
    """
    
    counts = getCounts(motifs)
    total = len(motifs) + 4
    for i in range(len(counts)):
        for nt in 'ACGT':
            counts[i][nt] += 1
            counts[i][nt] /= total
    return counts  # really now a profile
    
def getMostProbable(sequence, profile, k):
    """Return the most probable k-mer in sequence, given a profile."""
    
    maxProb = -1
    best_kmer = ''
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[j][kmer[j]]
        if prob > maxProb:
            maxProb = prob
            best_kmer = kmer
    return best_kmer
    
def getMostProbableMotifs(sequences, profile, k):
    """Return the most probable set of k-mers from a set of sequences,
       given a profile."""
       
    motifs = []
    for sequence in sequences:
        motifs.append(getMostProbable(sequence, profile, k))
    return motifs
    
def getConsensus(motifs):
    """Return the consensus sequence for a set of motifs."""
    
    counts = getCounts(motifs)
    
    consensus = ''
    for i in range(len(counts)):
        majority = 'A'
        for nt in counts[i]:
            if counts[i][nt] > counts[i][majority]:
                majority = nt
        consensus += majority
        
    return consensus
    
def getScore(motifs):
    """Return the score for a set of motifs."""
    
    counts = getCounts(motifs)
    consensus = getConsensus(motifs)
    
    t = len(motifs)
    score = 0
    for i in range(len(consensus)):
        nt = consensus[i]
        diff = t - counts[i][nt]
        score += diff
        
    return score



def readSequence(name):
    """
    read the sequences in a fasta file and return the list with each sequence as a string value
    Parameters
    ----------
    name : (string) the name of the FASTA file that contains the sequences

    Returns
    -------
    sequences: (list of strings)  list with each sequence as a string value 
    """

    record = SeqIO.index( (name + '.fasta'),'fasta')
    sequences = []
    for id in record:
        sequences.append(str(record[id].seq))

    return sequences                                   

def GibbsSampling(dna, k, t, N):
    motifs = []
    for i in range(t):
        index = random.randrange(len(dna[i])-k+1)
        motifs.append(dna[i][index:index+k])
    BestMotifs = motifs
    BestScore = getScore(motifs)
    for j in range(N):
        i = random.randrange(t)
        del motifs[i]
        profile = getProfileLaplace(motifs)
        prossibilities = {}
        sum = 0
        for x in range(len(dna[i])-k+1):
            candidate = dna[i][x:x+k]
            if candidate not in prossibilities:
                prossibilities[candidate] = 0
                for position in range(k):
                    nt = candidate[position]
                    prossibilities[candidate] += profile[position][nt]
            sum += prossibilities[candidate]
        pointer = 0
        possibility_keeper = 0
        chosen = random.random()
        chosen *= sum
        while (chosen > possibility_keeper):
            candidate = dna[i][pointer:pointer+k]
            possibility_keeper += prossibilities[candidate]
            pointer += 1
        motifs.insert(i,candidate)
        if getScore(motifs) < BestScore:
            BestMotifs = motifs
        return BestMotifs

def repeatGS(sequences,N):
    bestM = GibbsSampling(sequences, 20, len(sequences), 100000)
    bestS = getScore(bestM)
    for i in range(N-1):
        M = GibbsSampling(sequences, 20, len(sequences), 100000)
        S = getScore(M)
        if S<bestS :
            bestM = M
            bestS = S
    return(bestM)




def main():
    #section 2
    sequences = readSequence('upstream250')
   
    #section 1
    bestM = GibbsSampling(sequences, 20, len(sequences),10000000)
    bestS = getScore(bestM)
    print(bestM, bestS)

    #section 3
    bestM = repeatGS(sequences, 10000)
    bestS = getScore(bestM)
    print(bestM, bestS)
    
    #section 4 
    Consensus = getConsensus(bestM)
    print(Consensus)


main()