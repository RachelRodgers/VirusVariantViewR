# Stored in advance:
referenceCodon <- c("C", "A", "G")
names(referenceCodon) <- c(1,2,3)
referenceCodonString <- paste(referenceCodon, collapse = "")
originalAA <- Biostrings::GENETIC_CODE[[referenceCodonString]]

# Determined from the VCF file:
mutantAllele <- "C"
mutantPosition <- 3

# determine identity of mutant codon
mutantCodon <- referenceCodon
mutantCodon[mutantPosition] <- mutantAllele
mutantCodonString <- paste(mutantCodon, collapse = "")
mutantAA <- Biostrings::GENETIC_CODE[[mutantCodonString]]


## TEST ##


library("Biostrings")
View(GENETIC_CODE)
library("stringr")

setClass("Codon",
         representation = representation(sequence_vector = "character",
                                         protein = "character"),
         prototype = prototype(sequence_vector = NA_character_,
                               protein = NA_character_))

# read in the orf nucleotide sequence
orf3NT <- readr::read_file("../ORFs/Mod_CR6_ORF3_nt.txt")
orf3NTClean <- str_remove_all(orf3NT, pattern = "\\n")
#orf3NTVector <- unlist(base::strsplit(orf3NTClean, split = ""))

# split the orf nucleotide sequence into 3's, store in list
orf3Codons <- gsub("(.{3})", "-\\1", orf3NTClean)
orf3Codons <- stringr::str_split(orf3Codons, pattern = "-") #contains empty trailing "" that should be removed
orf3CodonsVec <- orf3Codons[[1]]
# remove empty element
orf3CodonsVec <- orf3CodonsVec[2:length(orf3CodonsVec)]
# put codons and start position in df
orf3CodonsDF <- data.frame("codon" = orf3CodonsVec, stringsAsFactors = FALSE)
orf3CodonsDF$position <- seq(from = 6681, to = 7307, by = 3)

# populate codon class
orf3CodonClassList <- vector(mode = "list", length = nrow(orf3CodonsDF))

for (i in 1:nrow(orf3CodonsDF)) {
  currentCodonInfo <- orf3CodonsDF[i, ]
  codonString <- unlist(base::strsplit(currentCodonInfo$codon, split = ""))
  startPosition <- currentCodonInfo$position
  names(codonString) <- seq(from = startPosition, to = startPosition + 2)
  
  currentCodonClass <- new("Codon", 
                           sequence_vector = codonString, 
                           protein = Biostrings::GENETIC_CODE[[currentCodonInfo$codon]])
  
  orf3CodonClassList[[i]] <- currentCodonClass
  names(orf3CodonClassList)[i] <- paste0("ORF3_codon_", i)
}

codonLUT <- rep(names(x = orf3CodonClassList), each = 3)
names(codonLUT) <- seq(from = 6681, to = 7307)

# How to find the correct codon to compare to?
mutantAllele <- "A"
mutantPosition <- 7250

referenceCodonName <- codonLUT[as.character(mutantPosition)]
# use the reference codon to determine what the mutated codon is:
referenceCodon <- orf3CodonClassList[[referenceCodonName]]

# determine identity of mutant codon
mutantCodon <- referenceCodon
mutantCodon@sequence_vector[as.character(mutantPosition)] <- mutantAllele
mutantCodonString <- paste(mutantCodon@sequence_vector, collapse = "")
mutantAA <- Biostrings::GENETIC_CODE[[mutantCodonString]]

# synonymous or non-synonymous ?
if (mutantAA != referenceCodon@protein) {
  print("non-synonymous")
} else {
  print("synonymous")
}

# if mutation > 1 character - print INDEL







codonInfo <- orf3CodonsDF[1, ]
# generate named sequence vector
codonString <- unlist(base::strsplit(codonInfo$codon, split = ""))
startPosition <- codonInfo$position
names(codonString) <- seq(from = startPosition, to = startPosition + 2)
codonString
currentProtein <- Biostrings::GENETIC_CODE[[codonInfo$codon]]

newCodon <- new("Codon",
                sequence_vector = codonString, protein = currentProtein)


# Get mutated codon information using "newCodon"
mutantAllele <- "C"
mutantPosition <- 6682














orf3CodonsList <- list(orf3CodonsVec)
sequenceIndex <- seq(from = 6681, to = 7307, by = 3)
names(orf3Codons) <- sequenceIndex

# populate the codon class







# name by position
# ORF3 positions
orf3Positions <- seq(from = 6681, to = 7307)
names(orf3NTVector) <- orf3Positions














namedORF3 <- orf3NTClean
names(orf3NTClean) <- orf3Positions

head(orf3NTClean)
tail(orf3NTClean)
