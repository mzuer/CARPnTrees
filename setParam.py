#######################################################################
################# Parameters used in the CARP program #################
# Marie Zufferey - UNIL - December 2015 ###############################
# License: Open Source AGPL v3 ########################################
#######################################################################

# You can edit this file to change some default parameters.

#************
# General parameters
#************
# needed to test if Windows or Python
from subprocess import Popen, PIPE
import platform

# Set working directory (needed for read and save files)
if platform.system() != "Windows":
    wd = Popen("pwd", shell=True, stdout=PIPE).communicate()[0].decode("ascii").strip()+"/"  # for UNIX
else:
    wd = ""

# Separator used in the input files (if only one input file, set the same for both)
sepP = ","                     # [char] Separator used in the polymorphisms file
sepH = ","                     # [char] Separator used in the haplogroups file

# Which column contains the haplogroup information
haploCol = 1                    # [integer] column that contains HG information (numbering not include row names)

# Which row to use as header ? (if only one input file, set the same for both)
hdrP = 0                    # [integer] row number(s) to use as the column names for the file with recent polymorphism
hdrH = 0                    # [integer] row number(s) to use as the column names for the file with haplogroups

# Format of the output
treePlotFormat = "png"          # [char] Format of the tree plot: png or svg
# in order to save tree 
treeFormat = "newick"           # [char] Format of the tree: Newick or NEXUS 
# ! Warning ! : if you change the treePlotFormat and treeFormat parameters,
# you need to comment/uncomment some lines in "CARPfcts.py" script
# (see explanations provided in this file)

plotName = "TreePlot"
treeName = "TreeText"

#************
# Parameters for Y-STR data only
#************
popCol = 0                      # [integer] In which column are the population name (#)
                                # set 0 if no column with population

polyCol1 = 2                    # [integer] first column that contains microsatellites data
polyCol2 = 0                    # [integer] last column that contains microsatellites data
                                # if microsatellite until from polyCol1 to the end of the file -> polyCol2 = 0

#************
# The following parameters are overwritten if the program is launched from commmand line
#************

# Verbose program or not ?
verb = False                                      # [boolean] True or False

# Comment one of the two following lines (...TreeP if paternal lineage, ...TreeM if maternal lineage)
# or give the path to your own file
allHaploTreeName = "allHaploTreeP.newick"        # for paternal lineage
# allHaploTreeName = "allHaploTreeM.newick"      # for maternal lineage

# Path to the files containing haplogroup and recent polymorphism data
if platform.system() != "Windows":
    polyFile = "foo.csv"             		 # [char] path to the file containing recent polymorphism data
    haploFile = "foo.csv"   			     # [char] path to the file containing the haplogroup data
else:
    polyFile = "data\\foo.csv"  	             # [char] path to the file containing recent polymorphism data
    haploFile = "data\\foo_HG.csv"               # [char] path to the file containing the haplogroup data

# Construct the tree using sequence data or microsatellite recent polymorphism data ?
analysis = "str"                # [char] "seq" or "str"
