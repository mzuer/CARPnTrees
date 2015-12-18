#######################################################################
################# Parameters used in the CARP program #################
# Marie Zufferey - UNIL - December 2015 ###############################
# License: Open Source AGPL v3 ########################################
#######################################################################
# You can edit this file to change some default parameters.

# Can be used e.g. to avoid to run the program from command line
# (in this case, comment some lines at the top of the script CARP.R, as explained in this latter file)

# By default "wd" is the current directory
# Can be changed e.g. to save outputed files in another directory
if(.Platform$OS.type == "unix"){
  wd <- system("pwd ", intern=TRUE)
} else{
  wd <- ""
}

# Name of the file that contains the whole haplogroup tree 
# (comment one of this two line or indicate the path to your file (must be in "newick" format)
allHaploTreeName <- "allHaploTreeM.newick"     # mtDNA haplogroups tree
# allHaploTreeName <- "allHaploTreeP.newick"   # Y-chromosome haplogroups tree

#************
# General parameters
#************

# Which column contains haplogroup information
haploCol <- 1  	                # [integer] column that contains HG (numbering not include row names)

# Separator used in the input files (if only one input file, set the same for both)
sepP <- ","                     # [char] Separator used in the recent polymorphism file
sepH <- ","                     # [char] Separator used in the haplogroups file

# Header in the files ? (if only one input file, set the same for both)
hdrP <- TRUE                    # [logical] Has the file with recent polymorphism data an header ?
hdrH <- TRUE                    # [logial] Has the file with haplogroups data an header ?

# Format of the output
treePlotFormat <- "png"         # [char] Format of the tree plot: png or svg
treeFormat <- "Newick"          # [char] Format of the tree: Newick or NEXUS
# ! Warning ! : if you change the treePlotFormat and treeFormat parameters,
# you need to comment/uncomment some lines in "CARPfcts.R" script
# (see explanations provided in this file)

# Name of the outputed tree
polyPlotname <- "TreePlot"	# [char] The name of your choice (used for saving picture of the tree)
polyTreename <- ""     		# [char] The name of your choice (used for saving textual representation of the tree)

#************
# If the program is run from command line, the following parameters are overwritten:
#************

# Name of the files containing haplogroup and recent polymorphism data
if(.Platform$OS.type == "unix"){
# if (Sys.info()["sysname"] != "Windows") 
  polyFile <- "data/myPolyData.csv"         
  haploFile <- "data/myHaploData.csv"
} else{
  polyFile <- "data\\myPolyData.csv"         # one of the "\" is the escape the "\" of the path name
  haploFile <- "data\\myHaploData.csv"  
}

# Construct the tree using mtDNA or Y-STR data ?
analysis <- "str"              # [char] "seq" or "str"

# By default, the program is not verbose
verb <- FALSE

