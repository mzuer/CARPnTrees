#! /usr/bin/python3.4
# -*- coding: utf-8 -*-

########################################################################
########## Script containing functions used in CARP programm ###########
# Marie Zufferey - UNIL - December 2015 ################################
# License: Open Source AGPL v3 #########################################
########################################################################
import sys
try:
    import numpy as np
    import pandas as pd
    import re
    import math
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    from Bio.Phylo.TreeConstruction import _DistanceMatrix
    from itertools import combinations
    from Bio import Phylo
    from optparse import OptionParser
#    import matplotlib.pyplot as plt # (only needed if saved as ".png")
except ImportError:
    sys.exit("Check you have the following packages correctly installed before running the script:\
             random, numpy, pandas, re, sys, math, itertools, optparse, dendropy,\
              Bio.Phylo, Bio.Align, Bio.AlignIO") # matplotlib.pyplot (if save as ".png")
from setParam import wd, sepP, haploCol, hdrP, treePlotFormat, treeFormat, analysis
from setParam import treeName, plotName
from CARPcls import *

markersList = list(np.genfromtxt("markersYSTR.csv", delimiter=",", usecols=1, skip_header=1, dtype=str))
ntList = list(np.genfromtxt("markersYSTR.csv", delimiter=",", usecols=3, skip_header=1, dtype=int))

#-------------------------------------------------------------------------------------------------------------#
# --------------------------- Function to retrieve information from command line -----------------------------#
#-------------------------------------------------------------------------------------------------------------#

def getOpt():
    '''
    Retrieve information and option from command line.
    :return:Retrieve information from command line (analysis="seq" or "str", PolyFile=file name, haploFile=file name)
    '''
    parser = OptionParser()
    parser.add_option("-f", "--format", action="store", type="string", dest="analysis", default="str",
                      help="type of analysis (str or seq)", )
    parser.add_option("-P", "--Poly", dest="polyFile", default=None, action="store", type="string",
                      help="complete name of polymorphism file (relative path from wd)")
    parser.add_option("-H", "--Haplo", dest="haploFile", default=None, action="store", type="string",
                      help="complete name of polymorphism file (relative path from wd), if not provided, are in Poly")
    parser.add_option("-v", "--verbose", dest="verb", default=False, action="store_true",
                      help="if omitted, nothing will be printed out in the terminal (default)")
    parser.add_option("-T", "--Tree", dest="tree", default=None, action="store", type="string",
                      help="name of the file with haplogroups or indicates if take default paternal or maternal tree")

    (options, args) = parser.parse_args()

    if options.tree == "m":
        options.tree = "allHaploTreeM.newick"      # the default maternal HG tree
    elif options.tree == "p":
        options.tree = "allHaploTreeP.newick"      # the default paternal HG tree
    elif options.tree == None:
        sys.exit("Error: missing \"-T\" argument.\nRead \"readme.txt\" could help.")
    elif not re.search(".newick$", options.tree):  # own HG tree can be provided in Newick format only
        sys.exit("Error: invalid \"-T\" argument.\nRead \"readme.txt\" could help.")

    if(options.analysis != "str" and options.analysis != "seq"):
        sys.exit("Error: invalid \"-f\" argument (\"str\" (default) or \"seq\".\nRead \"readme.txt\" could help.")

    if(not options.polyFile):
        sys.exit("Error: \"-P\" argument is mandatory.\nRead \"readme.txt\" could help.")

    if(not options.haploFile):
        options.haploFile = options.polyFile

    return (options.analysis, options.polyFile, options.haploFile, options.verb, options.tree)


#-------------------------------------------------------------------------------------------------------------#
# ------------------------------------------- Function to load data ------------------------------------------#
#-------------------------------------------------------------------------------------------------------------#


def loadData(haploFile, polyFile, analysis):
    '''
    Load the input data.
    :param haploFile: String: file name for haplogroups data.
    :param polyFile: String: file name for recent polymorphism data.
    :param analysis: String: type of data ("seq" for aligned sequences, "str" for microsatellites).
    :return: A GenData object (StrData or SeqData; class definition in CARPcls.py) representing the data for recent polymorphism.
    '''
    polyData=""
    if analysis == "str":
        data = pd.read_csv(wd + polyFile, index_col=0, header=hdrP, sep=sepP)
        if polyFile == haploFile:       # remove column with HG information
            data = data.drop(data.columns[[haploCol-1]], axis = 1)
        polyData = StrData(data)

    elif analysis == "seq":
        from Bio import AlignIO
        data = AlignIO.read(wd+polyFile, "fasta")
        polyData = SeqData(data)

    return polyData

#-------------------------------------------------------------------------------------------------------------#
# ------------------------------ Function to retrieve number of nts of microsat ------------------------------#
#-------------------------------------------------------------------------------------------------------------#


def ntMarker(marker):
    '''
    Retrieve the number of nucleotides of the repeat unit for a given Y-STR.
    :param marker: String: the name of the Y-STR.
    :return: Integer: the number of nucleotides of this marker, 1 if not known.
    '''

    # ntList and markersList defined in the other file
    marker = re.search("(^[A-Z]+[0-9]+)", marker).group(1)

    if marker in markersList:
        return ntList[markersList.index(marker)]
    else:
        return 1

#-------------------------------------------------------------------------------------------------------------#
# ----------------------------------- Functions to build and save the tree -----------------------------------#
#-------------------------------------------------------------------------------------------------------------#

####################
# Distance matrix
####################

def calcDist(data, labels, analysis):
    '''
    Compute the distance matrix for recent polymorphisms data.
    :param data: Recent polymorphism data (microsatellites data or aligned sequences)
                 for which we want calculate the distance matrix (markers=columns, individuals=rows).
    :param labels: List: the individuals ID (needed for row and column names of the distance matrix).
    :return: _DistanceMatrix object: the distance matrix (lower triangular part, without 0 diagonal).
    '''
    distM = ""

    if analysis == "str":
        markers = list(data)
        data = np.asmatrix(data)
        nts = [ntMarker(x) for x in markers]
        data = np.vstack([data, nts])    # add a row with nts value of microsat to the matrix
        allBruvo = []
        def calcBruvo (i1, i2, nts):
            diff = abs(i1-i2)/nts
            dist = 1-pow(2, -diff)
            return dist

        for ind1, ind2 in combinations(range(data.shape[0]-1), 2):
            bruvoD = 0

            for j in (range(data.shape[1])):
                nt = data[(data.shape[0]-1), j]   # nts is the last line of the matrix (number of nucleotide
                bruvoD += calcBruvo(data[ind1,j], data[ind2,j], nt)
            bruvoD /= (j+1)
            allBruvo.append(bruvoD)

        distM = _DistanceMatrix(labels)   # naming rows and columns with indiviudals ID

        index = 0

        for col in range(len(labels)-1):
            for row in range((1+col), len(labels)):
                distM[col, row] = allBruvo[index]
                index += 1

    elif analysis == "seq":
        from Bio.Phylo.TreeConstruction import DistanceCalculator
        from Bio.Align import MultipleSeqAlignment
        import itertools

        def pairwise(seq1, seq2):
            """
            Kimura 2-Parameter distance = -0.5*ln(1-2p-q) -0.25*ln(1-2q)
            p = percentage of transition
            q = percentage of transversion
            """
            from math import log, sqrt

            if len(seq1) != len(seq2):
                sys.exit("Not aligned sequences inputed")
            else:
                tot_len = len(seq1)

            transi = ["AG", "GA", "CT", "TC"]
            transv = ["AC","AT", "GC", "GT", "CA", "CG", "TA", "TG"]

            all_pairs = [i+j for i,j in zip(seq1, seq2)]

            transi_pairs = [x for x in all_pairs if x in transi]

            transv_pairs = [x for x in all_pairs if x in transv]

            p = float(len(transi_pairs)/tot_len)
            q = float(len(transv_pairs)/tot_len)


            a1 = 1-2*p-q
            a2 = 1-2*q

            try:
                k80_dist = -0.5*log(a1*sqrt(a2))  #k80_dist = -0.5*log(1-2*p-q) -0.25*log(1-2*q) # same 
            except ValueError:
                "error: cannot calculate log of negative values"

            return k80_dist


        def get_distance(msa):
           """Return a _DistanceMatrix for MSA object
           (taken from http://biopython.org/SRC/biopython/Bio/Phylo/TreeConstruction.py)
           :Parameters:
               msa : MultipleSeqAlignment
               DNA or Protein multiple sequence alignment.
           """
           if not isinstance(msa, MultipleSeqAlignment):
               raise TypeError("Must provide a MultipleSeqAlignment object.")

           names = [s.id for s in msa]
           dm = _DistanceMatrix(names)
           for seq1, seq2 in itertools.combinations(msa, 2):        # all possible combinations of pairs
               dm[seq1.id, seq2.id] = pairwise(seq1, seq2)
           return dm

        distM = get_distance(data)
        distM.names = labels

    return distM

#####################
# Build tree with NJ
#####################

def buildTree(distM):
    ''' Build the tree using NJ algorithm on distance matrix.
    :param distM: the distance matrix (_DistanceMatrix object).
    :return: Phylo.Tree build upon the distance matrix passed in argument using NJ algorithm.
    '''
    return DistanceTreeConstructor().nj(distM)


def buildTree_rooted_NJ(distance_matrix):
    """Construct and return an Neighbor Joining tree.
    :param distance_matrix: _DistanceMatrix object: the distance matrix on which tree is built.
    :return A Phylo.Tree object.

    Adapted from Yanbo Ye (2013), Bio.Phylo.DistanceTreeConstructor.nj()
    (only last part changed -> in this program, outputed nj tree must be rooted)

    This was the first version. Finally, the program uses an outgroup to root the tree.

    """

    if not isinstance(distance_matrix, _DistanceMatrix):
        raise TypeError("Must provide a _DistanceMatrix object.")

    from Bio.Phylo import BaseTree
    import copy

    # make a copy of the distance matrix to be used
    dm = copy.deepcopy(distance_matrix)
    # init terminal clades
    clades = [BaseTree.Clade(None, name) for name in dm.names]
    # init node distance
    node_dist = [0] * len(dm)
    # init minimum index
    min_i = 0
    min_j = 0
    inner_count = 0
    while len(dm) > 2:
        # calculate nodeDist
        for i in range(0, len(dm)):
            node_dist[i] = 0
            for j in range(0, len(dm)):
                node_dist[i] += dm[i, j]
            node_dist[i] = node_dist[i] / (len(dm) - 2)

        # find minimum distance pair
        min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
        min_i = 0
        min_j = 1
        for i in range(1, len(dm)):
            for j in range(0, i):
                temp = dm[i, j] - node_dist[i] - node_dist[j]
                if min_dist > temp:
                    min_dist = temp
                    min_i = i
                    min_j = j
        # create clade
        clade1 = clades[min_i]
        clade2 = clades[min_j]
        inner_count += 1
        inner_clade = BaseTree.Clade(None, "Inner" + str(inner_count))
        inner_clade.clades.append(clade1)
        inner_clade.clades.append(clade2)
        # assign branch length
        clade1.branch_length = (dm[min_i, min_j] + node_dist[min_i]
                                - node_dist[min_j]) / 2.0
        clade2.branch_length = dm[min_i, min_j] - clade1.branch_length


        # update node list
        clades[min_j] = inner_clade
        del clades[min_i]

        # rebuild distance matrix,
        # set the distances of new node at the index of min_j
        for k in range(0, len(dm)):
            if k != min_i and k != min_j:
                dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]
                                - dm[min_i, min_j]) / 2.0

        dm.names[min_j] = "Inner" + str(inner_count)
        del dm[min_i]

    # add a root and set last clade as one of the child and inner_clade (+ his descent) as another child
    root_name = "Inner" + str(inner_count+1)

    root = BaseTree.Clade(None, root_name)

    root.clades.append(clades[0])
    root.clades.append(clades[1])

    return BaseTree.Tree(root, rooted=True)


################
# Save the tree
################

def saveTree(tree, aux, verb=False):
    '''
    Save the tree in textual (Newick, NEXUS) or picture (png or svg) format.
    :param tree: A Phylo.Tree object: the tree for which to save the output.
    :param aux: String: an extra string that will end the name of the file (avoid overwritting during program execution).
    '''

    #import matplotlib.pyplot as plt (only needed if saved in ".png")

    # Save the tree in Newick format
    if treeFormat == "newick":
        ext = ".newick"

#    Phylo.write(tree, wd + treeName + aux + ext, treeFormat) # (*) uncomment this line if you want other than Newick format
# Written differently in order to get the tree in Newick format without any distance:
# If you want another format than Newick:
# - uncomment the above line marked with (*)
# - comment the following lines until the if verb statement marked with (*)
    filename = wd + treeName + aux + ext
    for i in tree.get_terminals():
        i.branch_length = 0
    for i in tree.get_nonterminals():
        i.branch_length = 0
    myTree = tree.format("newick")
    myTree = re.sub(":0.00000", "", myTree)

    orig_stdout = sys.stdout
    f = open(filename, 'w')
    sys.stdout = f
    print(myTree)
    sys.stdout = orig_stdout
    f.close()

    if verb:   # (*) see above: comment until here if you want another format than Newick
        print(treeFormat + " representation of the tree saved in \"" + wd + treeName + aux + ext + "\"")

    # Save plot of the tree (uncomment following lines to save tree as ".png")
    # you also need to uncomment the "import matplotlib.pyplot as plt" statement
    # (and this Python module must be installed on your computer)
    #Phylo.draw(tree, do_show=False)
    #plt.axis("off")
    #plt.savefig(wd + plotName + aux + "." + treePlotFormat)
    #if verb:
    #    print("Plot of the tree saved in \"" + wd + plotName + aux + "." + treePlotFormat + "\"")


def subsetPolyFile(samples, polyData):
    '''
    Subset from recent polymorphism data the data concerning the individuals passed in argument.
    :param samples: List of individuals for which to subset data.
    :param polyData: A GenData object containing the whole dataset from which it will be subset.
    :return: Recent polymorphism data concerning individuals only (if microsatellite data: data frame with STR data,
             if sequence data: AlignIO.MultipleSeqAlignment object)
    '''
    subPoly = ""
    if polyData.analysis == "str":
        subPoly = polyData.data.loc[samples]
    elif polyData.analysis == "seq":
        from Bio.Align import MultipleSeqAlignment
        pos = [polyData.labels.index(x) for x in samples]
        subPoly = MultipleSeqAlignment([])
        for i in range(len(pos)):
            subPoly.append(polyData.data[pos[i]])

    return subPoly

