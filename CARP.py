#! /usr/bin/python3.4
# -*- coding: utf-8 -*-

#######################################################################################
############################## Main script CARP program: ############################## 
############# building consensus tree from ancient and recent phylogenies #############  
# Marie Zufferey - UNIL - December 2015 ###############################################
# License: Open Source AGPL v3 ########################################################
#######################################################################################
import sys
try:
    from subprocess import Popen, PIPE
    import copy
    import dendropy
    import os
    import platform
    import random
    import itertools
    import math
    import pandas
    import re
    import optparse
#    import matplotlib.pyplot as plt # (uncomment if you want to save as the tree as ".png")
    from Bio import Phylo
    from Bio import AlignIO
    from Bio import Align
except ImportError:
    sys.exit("Check you have the following packages correctly installed before running the script:\
             random, numpy, pandas, re, sys, math, itertools, optparse, subprocess, dendropy,\
              Bio.Phylo, Bio.Align, Bio.AlignIO") # matplotlib.pyplot # if save as ".png"

from CARPfcts import *

# some "default" parameters are stored in another file
# so they can be modified easily
from setParam import wd, polyFile, haploFile, sepH, haploCol, hdrH, treeFormat, analysis, verb, allHaploTreeName

if platform.system() != "Windows":
    wd = Popen("pwd", shell=True, stdout=PIPE).communicate()[0].decode("ascii").strip()+"/"  # for UNIX
else:
    wd = ""

# Verbose or not (default: False - in setParam.py)
# (comment the following line if not run from command line)
verb = getOpt()[3]
if verb:
    print("#################################################################################")
    print("######### Building consensus tree for polymorphism and haplogroups data #########")
    print("#################################################################################\n")

# *********************************************************************************************************************
# ****************************** First build recent phylogenies using polymorphsim data *******************************
# *********************************************************************************************************************


###########################################
# Retrieve information from command line
###########################################

# if the script is run from command line, comment the next 4 lines
# and set parameters directly in setParam.py file
#
analysis = getOpt()[0]
polyFile = getOpt()[1]
haploFile = getOpt()[2]
allHaploTreeName = getOpt()[4]

if verb:
    print("******************** Parameters setting ********************")


inFormatPoly = re.search(".+\\.([a-z]+$)", polyFile).group(1)
inFormatHaplo = re.search(".+\\.([a-z]+$)", haploFile).group(1)

if analysis == "seq" and inFormatPoly != "fasta":
    sys.exit("Error: format for sequences must be FASTA file !\nRead \"readme.txt\" could help.")
elif analysis == "str" and inFormatPoly == "fasta":
    sys.exit("Errr: you did not choose the appropriate \"-f\" option for FASTA file.\nRead \"readme.txt\" could help.")
elif inFormatPoly == "fasta" and polyFile == haploFile:
    sys.exit("Error: if polymorphism data are in FASTA file, haplogroups data should be provided separately."+
             "\nRead \"readme.txt\" could help.")

##################
# Load the data
##################
if verb:
    print("\n... Data loading ...")
# Load the data for the polymorphisms (Y-STR and mtDNA HVR)
try:
    polyData = loadData(haploFile, polyFile, analysis)
except TypeError:
    sys.exit("Problem during data loading !\nRead \"readme.txt\" could help.")

if verb:
    print("Data succcessfully loaded.")
############################
# Compute distance matrix
############################
if verb:
    print("\n******************** Build recent polymorphism tree based on polymorphism data ********************")
    print("\n... Calculating matrix distance ...")
try:
    distM = calcDist(polyData.data, polyData.labels, polyData.analysis)

except:
    sys.exit("Problem during calculation of distance matrix !")
if verb:
    print("Distance matrix successfully computed.")


################################################
# Build the tree using NJ and save the result
################################################
if verb:
    print("\n... Building the tree using NJ algorithm ...")
try:
    polyTree = buildTree(distM)

except:
    sys.exit("Problem during building the tree from distance matrix !")

for i in polyTree.get_nonterminals():    # remove the "InnerX" node labels
    i.name = None

if verb:
    print("Tree successfully built.\n\n... Saving the results ...")
try:
    saveTree(polyTree, "Recent_py", verb)
except:
    sys.exit("Problem during results saving !")


# *********************************************************************************************************************
# **************************** Retrieve ancient phylogeny (HG) for individuals of the data ****************************
# *********************************************************************************************************************

##########
# small function to delete temporary file created in order to convert dendropy <-> Phylo Tree object
def rmTempTree():
    rm = "rm -f " + wd + "temp.newick"           # remove the "temp.newick" created above
    if platform.system() != "Windows":    # check if it works on Windows !
        os.system(rm)
##########
if verb:
    print("\n**************** From all haplogroups tree, select only haplogroups present in the data set ****************")

##################
# Load the data
##################
# Load the data for the haplogroups present in the dataset
try:
    haploData = pd.read_csv(wd + haploFile, index_col=0, header=hdrH, sep=sepH,
                            usecols=[0, haploCol])
except:
    sys.exit("Problem during loading the haplogroups data.\nSome parameters may need to be changed.\nRead \"readme.txt\" could help.")

haplo = list(haploData.iloc[:,0])  # list of HG present in the dataset

# Import the big tree with all HG
allHaploTree = dendropy.Tree.get(path=wd+allHaploTreeName, schema="newick")

###########################################################
# Nodes that contain individuals must be shifted as tips
###########################################################

# Get the nodes (not take the None)
allNodes = [x.label for x in allHaploTree.internal_nodes() if x.label]

# Nodes of the big haplogroups tree that represent haplogroups present in the data set
nodeLab = [x for x in allNodes if x in haplo]

# These nodes must be shifted as tips
for i in nodeLab:
    inode = allHaploTree.find_node_with_label(i)
    inode.label = None                           # replace label with None
    t = dendropy.Taxon(label=i)
    inode.new_child(taxon=t)                     # add branch with this label

for i in allHaploTree.leaf_node_iter():          # for all leaves label = taxon
    i.label = i.taxon

##############################################################################
# Remove from the big tree with all HG the ones that are not in the dataset
##############################################################################
if verb:
    print("\n... Subsetting from all haplogroups tree the part of interest ...")

allTips = [x.taxon.label for x in allHaploTree.leaf_node_iter() if x.taxon]

tipsToRemove = list(set(allTips)-set(haplo))

haploTree = allHaploTree
haploTree.prune_taxa_with_labels(tipsToRemove)          # work on a copy of allHaploTree

haploTree.write_to_path(wd + "temp.newick", treeFormat)   # save dendropy tree (for convert to Phylo)
haploTreeP = Phylo.read(wd + "temp.newick", treeFormat)

rmTempTree()

if verb:
    print("\n... Saving the haplogroups tree for the individuals present in the data set ...")
saveTree(haploTreeP, "Ancient_py", verb)

###############################
# Add samples ID to the tips
###############################

haploTreeID = copy.deepcopy(haploTree)

for x in haploTreeID.leaf_node_iter():
    individuals = haploData.loc[haploData.iloc[:,0] == x.taxon.label].index   # individuals sharing this haplogroup
    x.taxon.label += ":" + ",".join(individuals)

haploTreeID.write_to_path(wd + "temp.newick", treeFormat)
haploTreeIDP = Phylo.read(wd + "temp.newick", treeFormat)

rmTempTree()

if verb:
    print("\n... Saving the haplogroups tree with individuals ID as labels ...")
saveTree(haploTreeIDP, "AncientWithID_py", verb)

# *********************************************************************************************************************
# ********************** Build recent phylogeny and add it at the tip of the ancient phylogeny ************************
# *********************************************************************************************************************

######################################################
# For each haplogroup, subset data and build a tree
######################################################
if verb:
    print("\n******************** Build consensus tree joining polymorphism and haplogroups trees ********************")

uniqueHaplo = list(set(haplo))


for i in haploTree.leaf_nodes():
    hg = i.taxon.label

    samples = list(haploData[haploData.icol(0) == hg].index)

    if verb:
        print("\n... Building the tree for haplogroup "+hg+" shared by "+str(len(samples))+" individual(s) ...")

    if len(samples) == 1:
        lab = re.sub("\['", "", str(samples))
        lab = re.sub("'\]", "", lab)
        t = dendropy.Taxon(label=lab)
        i.new_child(taxon=t)                                  # add branch with indiviudals at this node

    elif len(samples) == 2:
        t1 = dendropy.Taxon(label=samples[0])
        t2 = dendropy.Taxon(label=samples[1])
        i.new_child(taxon=t1)
        i.new_child(taxon=t2)

    else:

        if i == haploTree.leaf_nodes()[0]:
            outHG = haploTree.leaf_nodes()[1].taxon.label
            outSample = list(haploData[haploData.icol(0) == outHG].index)
            nOut = random.randint(0, len(outSample)-1) #  randint(a,b) Return a random integer N such that a <= N <= b.
            outID = outSample[nOut]
        # add the outgroup to the subset
        samples.append(outID)
        subPoly = subsetPolyFile(samples, polyData)
        distM = calcDist(subPoly, samples, polyData.analysis)
        subPolyTree = buildTree(distM)
        Phylo.write(subPolyTree, "temp.newick", "newick")
        subPolyTree = dendropy.Tree.get(path="temp.newick", schema="newick")

                             # print(subPolyTree.as_ascii_plot())

        # To reroot a tree along an existing edge, you can use the reroot_at_edge method.
        # this method takes an edge object as as its first argument. this rerooting is a structural
        # change that will require the splits hashes to be updated before performing any tree comparisons
        # or calculating tree metrics. If needed, you can do this yourself by calling update_bipartitions
        node_out = subPolyTree.find_node_with_taxon_label(outID)
        subPolyTree.reroot_at_edge(node_out.edge, update_bipartitions=True)   # or False ???????????????????

        # set to None the taxon of the outgroup    # will always be single taxon as taken as outgroup ??????
        for x in subPolyTree.leaf_node_iter():
            if x.taxon.label == outID:
                x.taxon = None
                break
        subPolyTree.prune_leaves_without_taxa()
        i.add_child(subPolyTree.seed_node)
    # set the outgroup for next iteration:

    if len(samples) == 1:
        outID = samples[0]
    elif len(samples) == 2:
        nOut = random.randint(0, 1)
        outID = samples[nOut]
    else:
        # not take the last one because it is the outgroup
        nOut = random.randint(0, len(samples)-2)
        outID = samples[nOut]


for i in haploTree.leaf_nodes():
    i.edge_length = 0

rmTempTree()

haploTreeUltra = copy.deepcopy(haploTree)

for i in haploTree.internal_nodes():
    if isinstance(i.label, str) and re.search("Inner", i.label) != None:    #remove the default node labels "InnerX"
        i.label = None

haploTree.write_to_path(wd + "temp.newick", treeFormat)
haploTreeP = Phylo.read(wd + "temp.newick", treeFormat)

rmTempTree()

haploTreeUltra = copy.deepcopy(haploTree)

for i in haploTreeUltra.internal_nodes():
    i.label = None

for i in haploTreeUltra.internal_nodes():
    i.taxon = None
    i.edge_length  = None

orig_stdout = sys.stdout
f = open('ConsensusTreeWithoutNodLab_py.newick', 'w')
sys.stdout = f
print(a)
sys.stdout = orig_stdout
f.close()

if verb:
    print("\n... Saving the consensus tree...")
saveTree(haploTreeP, "ConsensusTree_py", verb)


if verb:
    if platform.system() == "Windows":
        print("\nDuring program execution, \"temp.newick\" was created. You can delete this file.")

    print("\n************************************ END ************************************")
    print("***************************** done with success *****************************")
    print("*****************************************************************************\n")
