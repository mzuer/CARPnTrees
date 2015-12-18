##### CARP PROGRAM #####
# with main function to launch the program with one function
########################
CARP_main_ancient <- function(polyFile, haploFile=NA, analysis="str", allHaploTreeName, 
                              headP=NA, sepP=NA, headH, sepH){
  
    source("CARPfcts.R")   # load functions used in this script
    source("CARPcls.R")
    wd=""
    verb=F
    haploCol = 1
    if(is.na(haploFile)){
        haploFile <- polyFile
    }
    
    loadPck("ape", verb)
    loadPck("seqinr", verb)
    
##### Subset from the whole haplogroups tree only the part with the haplogroups present in the dataset

    allHaploTree <- read.tree(paste0(wd,allHaploTreeName)) # load the tree with all HG
    haploTree <- allHaploTree
# Haplogroups present in the dataset
    haploData <- read.table(paste0(wd, haploFile), sep=sepH, header=headH, row.names=1)
    haplo <- unique(haploData[,haploCol])

# Nodes that contain individuals should be shifted as tips
    nodeLab <- haploTree$node.label[haploTree$node.label %in% haplo]
    for(i in nodeLab){
        tip <- list(edge=matrix(c(2,1),1,2), tip.label=i, edge.length=1.0, Nnode=1)
        class(tip) <- "phylo"
        haploTree$node.label[haploTree$node.label==i] <- "temp"
        haploTree <- bind.tree(haploTree, tip, where=findPos("temp", haploTree))
        haploTree$node.label[haploTree$node.label=="temp"] <- ""
    }

##### Select from the tree only the branches with the haplogroups of the data
# Tips of the big tree that are not present in the dataset, remove it
            # Not remove internal nodes because they could contain individuals
            # Repeat until all tips contain individuals
    allTips <- haploTree$tip.label
    tipsToRemove <- setdiff(allTips, haplo)

    haploTree <- drop.tip(haploTree, tipsToRemove, trim.internal = T)

##### Add individuals ID at the tips and nodes of the tree
# work on a copy of haploTree
    haploTreeID <- haploTree

    addID <- function(HG){   
  # Small function that takes the haplogroup name in argument and returns a string
  # with the collapsed ID of the individuals that belong to this haplogroup
        ID <- paste0(rownames(haploData)[which(haploData[,haploCol]==HG)], collapse="\n")
        return(ID)
    }

# Add ID to the tips
    tipsHG <- haploTreeID$tip.label     # tips with individuals from the dataset
    tipsID <- sapply(tipsHG, FUN=addID)
haploTreeID$tip.label <- paste(tipsHG, tipsID, sep="_")

  # haploTreeID$tip.label <- mixedFontLabel(tipsHG, tipsID, bold = 1, italic = 2, sep=" ")

  haploTreeID$edge.length <- rep(0, length(haploTreeID$edge.length))
  newickTree <- write.tree(haploTreeID)
  newickTree <- gsub(":0", "", newickTree)

    return (newickTree)

}
#############################################################################################################
CARP_plot_ancient <- function(polyFile, haploFile=NA, analysis="str", allHaploTreeName, 
                              headP=NA, sepP=NA, headH, sepH){
  
  source("CARPfcts.R")   # load functions used in this script
  source("CARPcls.R")
  wd=""
  verb=F
  haploCol = 1
  if(is.na(haploFile)){
    haploFile <- polyFile
  }
  
  loadPck("ape", verb)
  loadPck("seqinr", verb)
  
  
  ##### Subset from the whole haplogroups tree only the part with the haplogroups present in the dataset
  
  allHaploTree <- read.tree(paste0(wd,allHaploTreeName)) # load the tree with all HG
  haploTree <- allHaploTree
  # Haplogroups present in the dataset
  haploData <- read.table(paste0(wd, haploFile), sep=sepH, header=headH, row.names=1)
  haplo <- unique(haploData[,haploCol])
  
  # Nodes that contain individuals should be shifted as tips
  nodeLab <- haploTree$node.label[haploTree$node.label %in% haplo]
  for(i in nodeLab){
    tip <- list(edge=matrix(c(2,1),1,2), tip.label=i, edge.length=1.0, Nnode=1)
    class(tip) <- "phylo"
    haploTree$node.label[haploTree$node.label==i] <- "temp"
    haploTree <- bind.tree(haploTree, tip, where=findPos("temp", haploTree))
    haploTree$node.label[haploTree$node.label=="temp"] <- ""
  }
  
  ##### Select from the tree only the branches with the haplogroups of the data
  # Tips of the big tree that are not present in the dataset, remove it
  # Not remove internal nodes because they could contain individuals
  # Repeat until all tips contain individuals
  allTips <- haploTree$tip.label
  tipsToRemove <- setdiff(allTips, haplo)
  
  haploTree <- drop.tip(haploTree, tipsToRemove, trim.internal = T)
  
  ##### Add individuals ID at the tips and nodes of the tree
  # work on a copy of haploTree
  haploTreeID <- haploTree
  
  addID <- function(HG){   
    # Small function that takes the haplogroup name in argument and returns a string
    # with the collapsed ID of the individuals that belong to this haplogroup
    ID <- paste0(rownames(haploData)[which(haploData[,haploCol]==HG)], collapse="\n")
    return(ID)
  }
  
  # Add ID to the tips
  tipsHG <- haploTreeID$tip.label     # tips with individuals from the dataset
  tipsID <- sapply(tipsHG, FUN=addID)
#   haploTreeID$tip.label <- paste(tipsHG, tipsID, sep="_")
  
  haploTreeID$tip.label <- mixedFontLabel(tipsHG, tipsID, bold = 1, italic = 2, sep=" ")
  
  plot.phylo(haploTreeID)

#   haploTreeID$edge.length <- rep(0, length(haploTreeID$edge.length))
#   newickTree <- write.tree(haploTreeID)
#   newickTree <- gsub(":0", "", newickTree)
#   
#   return (newickTree)
  
}
