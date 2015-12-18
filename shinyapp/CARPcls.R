######################################################################
# Classes used in CARP scripts                                       #
######################################################################

loading <- try(library("methods", verbose=F))
if(class(loading)=="try-error"){
  loading2 <- try(install.packages("methods", verbose=F))
    if(class(loading2)=="try-error"){
      stop("\"methods\" package is needed and could not be installed.")      
    }
}  


loading <- try(library("seqinr", verbose=F))
if(class(loading)=="try-error"){
  loading2 <- try(install.packages("methods", verbose=F))
    if(class(loading2)=="try-error"){
      stop("\"methods\" package is needed and could not be installed.")      
    }
}

setClass("genData", 
         representation(
           labels = "character",
           type = "character")
)

setMethod("show", "genData",
          function(x){
            b <- ifelse(x@type=="str", "STR", "sequence")
            cat(paste0("Contains ", b, " data used to build recent phylogeny.\n"))
            cat(paste0("Individuals: "));cat(x@labels)
            if(x@type=="str"){
              cat(paste0("\nMarkers: "));cat(x@markers)
            }
            cat("\n")
          }
)

setClass("strData", 
         representation(
           data =  "data.frame",
           markers = "character"),
         contains="genData"
)

setClass("seqData", 
         representation(
           data =  'ANY'),
         contains="genData"
)
