# First, we want to look at the number of mutations at a particular residue site
# Second, we want to look at the porportion of mutations that don't lead to a mutation per residue site 
# Third, we want to look at the number of different barcodes per mutation (for the same site, and then again for codon change.)
# Fourth, we want to look at most favored amino acid substitution across the whole protein sequence

import pickle 
import numpy as np
import matplotlib.pyplot as plt
import seaborn 
import matplotlib
import copy

bar2codon = pickle.load(open("allele_dic.pkl","rb"))
codon2aa = pickle.load(open("translate.pkl","rb"))
aa2num = pickle.load(open("aminotonumber.pkl", "rb"))


def mutationspersite():
    mutationspersite = np.zeros(77)
    datamat = np.zeros((21, 77));
    #i have mutations per site
    #barcodes to amino acid
    AAs = {};
    for key in bar2codon:
        tmp = bar2codon[key][0].split("_") #locate mutation site
        #take codon, map which amino acid
        tmp[1] = tmp[1].replace('T', 'U')
        #print(tmp[1])
        datamat[ aa2num[codon2aa[tmp[1]]]][float(tmp[0])-1] += 1
	mutationspersite[float(tmp[0])-1] += 1 #tally in site
    x_axis = np.arange(1,78)
    #plt.scatter(x_axis, mutationspersite)
    datamatnozero = copy.deepcopy(datamat)
    datamatnozero[datamatnozero == 0] = 0.7
    datamatlog = np.log10(datamatnozero);

    
    heatmap = plt.pcolor(datamatlog, cmap='PiYG')
    #print(datamat[10][10])
    #print(np.log(datamat[10][10]))
    plt.colorbar(heatmap);
    print(np.min(datamatlog), np.max(datamatlog))
    matplotlib.colors.Normalize(vmin=np.min(datamatlog), vmax=np.max(datamatlog), clip=True)
    matplotlib.colors.Normalize(vmin=-1, vmax=np.max(datamatlog), clip=True)
    #matplotlib.colors.Normalize.autoscale(datamatlog)
    plt.title("Heatmap")
    plt.ylabel("Amino Acids")
    plt.xlabel("Residue Number")
    plt.show()

def barcodes_per_AA():
    #get number of barcodes that map to the same codon
    #i might need number o funique codons
    #unique_codon = {};
    #per residue, number of barcodes that map to 
    AAs = {};
    for key in bar2codon:
        tmp = bar2codon[key][0].split("_") #locate mutation site
        AAs[ codon2aa[tmp[1]] ] = 0;
    for key in bar2codon:
        tmp = bar2codon[key][0].split("_") #locate mutation site
        AAs[ codon2aa[tmp[1]] ] += 1;
        
    stuff = plt.scatter(np.arange(1, 21), AAs)
    plt.show();






def main():
    mutationspersite()
#    barcodes_per_AA();


if __name__ == "__main__":
    main(); 
