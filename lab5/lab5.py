"""
CS 309 - Computational Biology
Lab 5
Brian Nguyen, Kairuo Yan
"""

from re import T
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt

"""
This function reads the Fasta file, coverts the content into a Python dictionary with the keys being the headers and 
the values being the genomes

parameters: filename: name of the Fasta file
ptype string
return: sequenceDict: a Python dictionary with the keys being the headers and the values being the genomes
rtype dictionary
"""
def getDict(filename):
    sequenceDict = {rec.id : str(rec.seq) for rec in SeqIO.parse(filename, "fasta")}
    return sequenceDict


"""
This function finds the genomes that have the B.1.1.7 variant in one sample and their proportion in the said sample

parameters: sequenceDict: a Python dictionary with the keys being the headers and the values being the genomes
            is0512: a boolean value checking if the sample is from the date 05 of december 2020
ptype dictionary
      boolean
return: proportion: the proportion of genomes that have the B.1.1.7 variant in the sample
        labels: the list containing the genomes that have the B.1.1.7 variant in the sample 
rtype double
      list
"""
def findProportion(sequenceDict, is0512):
    count = 0
    labels = []
    N501pos = 23063
    H69_70pos = 21765
    P681pos = 23604
    if is0512:
        N501pos = 23071
        H69_70pos = 21764
        P681pos = 23612
    for key in sequenceDict:
        if sequenceDict[key][N501pos-1]=='T' and sequenceDict[key][H69_70pos-1:H69_70pos+5]=='------' and sequenceDict[key][P681pos-1]=='A':
            count += 1
            labels.append(key)
    proportion = (count / (len(sequenceDict)-1)*100)
    return proportion, labels


"""
This function makes the bar graph showing the proportions of genomes having the B.1.1.7 variants for all 
five sampling dates and returns all the genomes that have the B.1.1.7 variants for all five sampling dates

parameters: fastaFiles: the list containing the names of sampling date files
ptype list
return: genomesWithB117Variant: the list containing all the genomes that have the B.1.1.7 variants for all five sampling dates
rtype list
"""
def makeBarGraph(fastaFiles):
    proportions = []
    genomesWithB117Variant = []
    for file in fastaFiles:
        is0512 = False
        if '0512' in file:
            is0512 = True
        sequenceDict = getDict(file)
        proportion, labels = findProportion(sequenceDict, is0512)
        proportions.append(proportion)
        genomesWithB117Variant.append(labels)
    samples = [i[0:6] for i in fastaFiles]
    y_pos = np.arange(len(samples))
    plt.bar(y_pos, proportions, align = 'center', alpha = 0.5)
    plt.xticks(y_pos, samples)
    plt.ylabel('Proportion as percentage')
    plt.title('Proportions of genomes with B.1.1.7 variants for each sample date')
    plt.show()
    return genomesWithB117Variant


"""
This function reads the .pim file, coverts the content of each line into a Python string, puts all said lines into a list, and 
returns the list

parameters: fileName: name of the .pim file
ptype string
return: lineStrings: the list containing the content of the file
rtype list
"""
def get_file_lines(fileName):
    lineStrings = []
    FILE = open(fileName, "r")
    for line in FILE:
        if (line[-1] == "\n"): lineStrings.append(line[:-1])
        else: lineStrings.append(line)
    FILE.close()
    return lineStrings


"""
This function converts a list containing the content of the .pim file, converts it into a python dataframe, and returns the dataframe 
with names of the genomes stored in a list 

parameters: lines: the list containing the content of the file
ptype list
return: DF: the dataframe representing the file
        names: the list containing the names of the genomes in the file
rtype DataFrame
      list
"""
def lines_to_dataframe(lineStrings):
    names = []
    dist_dict = dict()
    for line in lineStrings:

        linarr = line.split()
        percents = list(map(lambda x : float(x), linarr[2:]))
        anno = linarr[1].split("|")[1]

        if anno in dist_dict.keys(): print(anno)
        dist_dict[anno] = percents
        names.append(anno)

    DF = pd.DataFrame(columns = names)
    for sample in names: DF.loc[sample] = dist_dict[sample]
    return DF, names


"""
This function converts a .pim file into a python dataframe by calling the 2 helper functions above

return: DF: the dataframe representing the file
        names: the list containing the names of the genomes in the file
rtype DataFrame
      list
"""
def pim_to_dataframe(fileName : str) -> pd.DataFrame:
    lines = get_file_lines(fileName)
    return lines_to_dataframe(lines)


"""
This function creates a distance matrix and a list of taxa from a .pim file

parameters: filename: name of the .pim file
ptype string
return: distanceMatrix: the 2d array represents the distance matrix 
        taxa: the list containing the taxa
rtype nparray
      list
"""
def makeDistanceMatrixAndTaxaList(fileName):
    df, taxa = pim_to_dataframe(fileName)
    distanceMatrix = df.to_numpy()
    for r in range(len(distanceMatrix)):
        for c in range(len(distanceMatrix)):
            distanceMatrix[r][c] = 100 - distanceMatrix[r][c]

    return distanceMatrix, taxa


"""
This function computes all the distances from all the taxa to another and contains
the results in the list

parameters: taxa: the list containing the taxa
            distanceMatrix: the 2d array represents the distance matrix
ptype list
      nparray
return: R: the list containing all the distances from all the taxa to another
rtype nparray
"""
def computeR(distanceMatrix, n):
    R = np.zeros(n)
    for i in range(n):
        for j in range(n):
            R[i] += distanceMatrix[i,j]
    return R


"""
This function makes matrix M 

parameters: R: the list containing all the distances from all the taxa to another
            distanceMatrix: the 2d array represents the distance matrix
            taxa: the list containing the taxa
ptype nparray
      nparray
      list
return: M: the matrix M
rtype nparray
"""
def makeMatrixM(distanceMatrix, R, n):
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i != j:
                M[i,j] = (n-2) * distanceMatrix[i,j] - R[i] - R[j]
    return M


"""
This function makes matrix M 

parameters: R: the list containing all the distances from all the taxa to another
            distanceMatrix: the 2d array represents the distance matrix
            taxa: the list containing the taxa
ptype nparray
      nparray
      list
return: M: the matrix M
rtype nparray
"""
def modifyDistanceMatrix(distanceMatrix, i, j, taxa, R):
    elementI = taxa[i]
    elementJ = taxa[j]
    Dv = np.zeros(len(distanceMatrix)+1)   #list containing the distance between the new node V and the rest of the taxa in the distance matrix
    for k in range(len(taxa)):
        if (k!=i) and (k!=j):
            Dvk = (distanceMatrix[i,k] + distanceMatrix[j,k] - distanceMatrix[i,j]) /2    #the distance between the new node v and k
            Dv[k] = Dvk
    
    Lvi = ( (distanceMatrix[i,j] / 2) + ((R[i]-R[j])/ (2* (len(distanceMatrix)-2) )) )
    Lvj = ( (distanceMatrix[i,j] / 2) + ((R[j]-R[i])/ (2* (len(distanceMatrix)-2) )) )

    Vnode = '('+ elementI+':'+ str(Lvi) +','+ elementJ  + ':'+ str(Lvj) +')'

    # np.insert(len(distanceMatrix), Vnode, Dv[:-1])

    # D.insert(len(D), V_name, D_v[:-1])
    # D.loc[V_name] = D_v

    distanceMatrix = np.delete(distanceMatrix,i,0)
    distanceMatrix = np.delete(distanceMatrix,j,0)
    distanceMatrix = np.delete(distanceMatrix,i,1)
    distanceMatrix = np.delete(distanceMatrix,j,1)



    return distanceMatrix, elementI, elementJ, Vnode


"""
This function performs the neighbor joining algorithm

parameters: distanceMatrix: the 2d array represents the distance matrix 
        taxa: the list containing the taxa
rtype nparray
      list
return: newickString: the newick string represents the phylogenetic tree 
rtype string
"""
def neighborJoining(distanceMatrix, taxa):
    while(len(distanceMatrix))>2:
        n = len(taxa)
        R = computeR(distanceMatrix, n)
        M = makeMatrixM(distanceMatrix, R, n)
        minFlatIndex = np.argmin(M)
        (i,j) = np.unravel_index(minFlatIndex, (n,n))
        distanceMatrix, elementI, elementJ, Vnode = modifyDistanceMatrix(distanceMatrix, i, j, taxa, R)
        taxa.remove(elementI)
        taxa.remove(elementJ)
        taxa.append(Vnode)
    newickString = '('+ taxa[0] +',' + taxa[1] + ':' + str(distanceMatrix[0,1]) +')'
    return newickString


"""
This function visualzes the newick string in a tree format

parameters: distanceMatrix: the 2d array represents the distance matrix 
        taxa: the list containing the taxa
rtype nparray
      list
return: none
"""
# def visualizeTree(distanceMatrix, taxa):
#     newickString = neighborJoining(distanceMatrix, taxa)
#     tree = Phylo.read(StringIO(newickString), 'newick')
#     Phylo.draw(tree)


def main():
    # fastaFile = 'HCOV19-ENGLAND-081220-A.fasta'
    fastaFile = '081220_fasta_alignment.fasta'
    pimFile = 'HCOV19-ENGLAND-051220-D.pim'
    fastaFiles = ['031120_fasta_alignment.fasta', '101120_fasta_alignment.fasta', '171120_fasta_alignment.fasta', '051220_fasta_alignment.fasta', '081220_fasta_alignment.fasta']
    # genomesWithB117Variant = makeBarGraph(fastaFiles)
    # print(genomesWithB117Variant)
    distanceMatrix, taxa = makeDistanceMatrixAndTaxaList(pimFile)
    str = neighborJoining(distanceMatrix, taxa)
    # print(distanceMatrix)
    # print(len(taxa))
    # R = computeR(distanceMatrix, taxa)
    # print(makeMatrixM(distanceMatrix, R, taxa))
    #a=makeDistanceMatrixAndTaxaList(pimFile)[1]
    # print(len(distanceMatrix))
    # visualizeTree(distanceMatrix, taxa)
    #print(a)
    
main()