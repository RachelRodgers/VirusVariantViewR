# Mod_CR6_ORF_Information.R

# Generate codon classes and look-up-tables for each ORF in ModCR6.

#----- Load Libraries -----#
library("Biostrings")
library("stringr")

#----- Defininitions -----#
setClass("Codon",
         representation = representation(sequence_vector = "character",
                                         protein = "character"),
         prototype = prototype(sequence_vector = NA_character_,
                               protein = NA_character_))

PopulateCodonClasses <- function(orfFile, orfName, orfStart, orfEnd) {
  # read in, clean & separate ORFs
  orfNTRaw <- readr::read_file(orfFile)
  orfNTNoBreaks <- str_remove_all(orfNTRaw, pattern = "\\n")
  orfCodonsList <- stringr::str_split(gsub("(.{3})", "-\\1", orfNTNoBreaks),
                                      pattern = "-")
  orfCodonsVec <- orfCodonsList[[1]]
  orfCodonsVec <- orfCodonsVec[2:length(orfCodonsVec)] # remove empty element
  
  # put codons and start positions in df
  orfCodonsDF <- data.frame("codon" = orfCodonsVec, stringsAsFactors = FALSE)
  orfCodonsDF$position <- seq(from = orfStart, to = orfEnd, by = 3)
  
  # populate codon class list
  orfCodonClassList <- vector(mode = "list", length = nrow(orfCodonsDF))
  
  for(i in 1:nrow(orfCodonsDF)) {
    currentCodonInfo <- orfCodonsDF[i, ]
    codonString <- unlist(base::strsplit(currentCodonInfo$codon, split = ""))
    currentStartPos <- currentCodonInfo$position
    names(codonString) <- seq(from = currentStartPos, 
                              to = currentStartPos + 2) # because codons come in 3's
    currentCodonClass <- new("Codon",
                             sequence_vector = codonString,
                             protein = Biostrings::GENETIC_CODE[[currentCodonInfo$codon]])
    orfCodonClassList[[i]] <- currentCodonClass
    names(orfCodonClassList)[i] <- paste0(orfName, "_codon_", i)
  }
  return(orfCodonClassList)
}

#----- Generate ORF Codon Classes -----#

# ~ ORF1 ~ #
# Positions: 6 - 5069
# Length: 5064 nt | 1687 aa
orf1CodonClassList <- PopulateCodonClasses(orfFile = "../Mod_CR6_ORFs/Mod_CR6_ORF1_nt.txt",
                                           orfName = "ORF1",
                                           orfStart = 6, orfEnd = 5069)
orf1CodonLUT <- rep(names(x = orf1CodonClassList), each = 3)
names(orf1CodonLUT) <- seq(from = 6, to = 5069)

# ~ ORF2 ~ #
# Positions 5056 - 6681
# Length: 1626 nt | 541 aa
orf2CodonClassList <- PopulateCodonClasses(orfFile = "../Mod_CR6_ORFs/Mod_CR6_ORF2_nt.txt",
                                           orfName = "ORF2",
                                           orfStart = 5056, orfEnd = 6681)
orf2CodonLUT <- rep(names(x = orf2CodonClassList), each = 3)
names(orf2CodonLUT) <- seq(from = 5056, to = 6681)

# ~ ORF3 ~ #
# Positions: 6681 - 7307
# Length: 627 nt | 208 aa
orf3CodonClassList <- PopulateCodonClasses(orfFile = "../Mod_CR6_ORFs/Mod_CR6_ORF3_nt.txt",
                                           orfName = "ORF3",
                                           orfStart = 6681, orfEnd = 7307)

orf3CodonLUT <- rep(names(x = orf3CodonClassList), each = 3)
names(orf3CodonLUT) <- seq(from = 6681, to = 7307)

# Biostrings AMINO_ACID_CODE with Stop Codon Added
aminoAcidCode <- c(Biostrings::AMINO_ACID_CODE, "*" = "*")


save.image("Mod_CR6_ORF_Information.RData")
