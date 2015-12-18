#################################################################
################### CARPnTrees program readme ###################
#################################################################

************* Before running the CARP script *************

Check the 1st line of this script: shebang may need to be changed.

By default, program should be launched from the directory containing the R and python scripts.

************* Arguments of the command line *************

The programm is launched in a terminal as followed (.* = .py or .R):

./CARP.* [-f, -v] -P FILE_POLY_HAPLO -T FILE_HG_TREE            

./CARP.* [-f, -v] -H FILE_HAPLO -P FILE_POLY -T FILE_HG_TREE

Mandatory argument:
- the file containing the recent polymorphism data following the "-P" or "--Poly" flag
- indication about the file to use for the haplogroups tree after the "-T" or "--Tree" flag, it can be:
      -T m         => for maternal lineages, use the default haplogroups tree provided with the program
      -T p   	   => for paternal lineages, use the default haplogroups tree provided with the program [not yet available !]
      -T *.newick  => use the haplogroups tree specified after the flag (must be Newick format, ".newick")

Optional arguments:
-H, --Haplo: the file containing the haplogroups data following the "-H" flag.
If no file is provided, haplogroups and polymorphisms data are in the same file.

-f, --format: must be (without quotes) "str" (default) or "seq". Indicates if the polymorphism data
are aligned FASTA sequences ("seq") or STR data ("str").

-v, --verbose: if omitted, nothing will be printed out in the terminal (default).

If you don't want to run the script from command line, set some parameters in the setParam.* file
and comment ~6 lines in CARP.* file (where the function "getOpt" is called - see beginning of CARP.*).


************* Shape of the data *************

By default, the data must be formated as follows:
- first column (row names) must be sample ID
- haplogroups are in the 2nd column (1st column after row names)
  and sample ID the first one (row names)
- if polymorphisms data are aligned sequences, they must be in FASTA format
- sequences from FASTA file must be aligned (i.e. same length)
- the first line of the files (except FASTA) contains an header
- default separator is ","

Most of these parameters and other ones (e.g. name of output files, or if your data are not 
in the same directory as the script) can easily be changed in the file "setParam.*". 


************* Version and requirements *************

If needed R packages are not installed, they are installed during the program execution.
R script tested with versions: 
- R 3.2.2
- ape 3.3
- seqinr 3.1.4
- optparse 1.3.2
- methods 3.2.2

Python required modules should be installed before program launching.
Python script tested with versions: 
- Bio.Phylo (Biopython 1.66)
- dendropy 4.0.3
- Bio.AlignIO (Biopython 1.66)
- Bio.Align (Biopython 1.66)
- numpy 1.8.2
- pandas 0.13.1
- matplotlib 1.3.1
- sys 3.4.3
- re 2.2.1
- itertools
- optparse
- subprocess
- copy
- os
