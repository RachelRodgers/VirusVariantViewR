# CW3_ORF_Information.R

# Generate codon classes and look-up-tables for each ORF in ModCR6.

#----- Load Libraries -----#
library("Biostrings")
library("tidyverse")

#----- Defininitions -----#
setClass("Codon",
         representation = representation(sequence_vector = "character",
                                         protein = "character",
                                         gene_name = "character",
                                         protein_position = "character"),
         prototype = prototype(sequence_vector = NA_character_,
                               protein = NA_character_,
                               gene_name = NA_character_,
                               protein_position = NA_character_)
         )

setClass("ORF",
         representation = representation(name = "character",
                                         range = "character",
                                         LUT = "character",
                                         codons = "list"),
         prototype = prototype(name = NA_character_,
                               range = NA_character_,
                               LUT = NA_character_))

PopulateCodonClasses <- function(orfFile, orfName, orfStart, orfEnd) {
  
  # data needed for numbering the AA position of ORF1 proteins:
  # add protein name and position - varies for ORF1 by position, all others constant
  # Set up function to "re-set" amino acid position in genes in ORF1.  We want
  #   to know amino-acid positions with couting re-starting at the first aa in
  #   the protein, not the consecutive nt number or consecutive aa number (i)
  #   except for the first gene, NS1/2.
  orf1Genes <- c("NS3", "NS4", "NS5", "NS6", "NS7")
  ResetAANum <- function(startNT, stopNT) {
    # start - starting nt position of gene
    # end -ending nt position of gene
    # generate sequence 1...n
    seqVector <- as.character(seq_along(startNT:stopNT))
    names(seqVector) <- startNT:stopNT
    return(seqVector)
  }
  orf1GenesList <- pmap(.l = list(startNT = list("NS3" = 342,
                                                 "NS4" = 706,
                                                 "NS5" = 871,
                                                 "NS6" = 995,
                                                 "NS7" = 1178),
                                  stopNT = list("NS3" = 705,
                                                "NS4" = 870,
                                                "NS5" = 994,
                                                "NS6" = 1177,
                                                "NS7" = 1688)),
                        .f = ResetAANum)
  
  # read in, clean & separate ORFs
  orfNTRaw <- readr::read_file(orfFile)
  orfNTNoBreaks <- str_remove_all(orfNTRaw, pattern = "[:cntrl:]")
  orfCodonsList <- stringr::str_split(gsub("(.{3})", "-\\1", orfNTNoBreaks),
                                      pattern = "-") # split into groups of 3 letters
  orfCodonsVec <- orfCodonsList[[1]] # put groups in vector
  orfCodonsVec <- orfCodonsVec[2:length(orfCodonsVec)] # remove empty element @ beginning of vec
  
  # put codons and start positions in df
  orfCodonsDF <- data.frame("codon" = orfCodonsVec, stringsAsFactors = FALSE) # start as single-row df
  orfCodonsDF$position <- seq(from = orfStart, to = orfEnd, by = 3) # add numbers to df
  
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
    
    if (orfName == "ORF1") {
      currentCodonClass@gene_name <- case_when(i %in% 1:341 ~ "NS1/2",
                                               i %in% 342:705 ~ "NS3",
                                               i %in% 706:870 ~ "NS4",
                                               i %in% 871:994 ~ "NS5",
                                               i %in% 995:1177 ~ "NS6",
                                               i %in% 1178:1688 ~ "NS7")
      # Re-set the start of the protein so numbering begins at gene start.
      currentCodonClass@protein_position <- case_when(i %in% 1:341 ~ as.character(i),
                                                      i %in% 342:705 ~ orf1GenesList$NS3[as.character(i)],
                                                      i %in% 706:870 ~ orf1GenesList$NS4[as.character(i)],
                                                      i %in% 871:994 ~ orf1GenesList$NS5[as.character(i)],
                                                      i %in% 995:1177 ~ orf1GenesList$NS6[as.character(i)],
                                                      i %in% 1178:1688 ~ orf1GenesList$NS7[as.character(i)])
    } else if (orfName == "ORF2") {
      currentCodonClass@gene_name <- "VP1"
      currentCodonClass@protein_position <- as.character(i)
    } else if (orfName == "ORF3") {
      currentCodonClass@gene_name <- "VP2"
      currentCodonClass@protein_position <- as.character(i)
    } else if (orfName == "ORF4") {
      currentCodonClass@gene_name <- "VF1"
      currentCodonClass@protein_position <- as.character(i)
    }
    
    orfCodonClassList[[i]] <- currentCodonClass
    names(orfCodonClassList)[i] <- paste0(orfName, "_codon_", i)
  }
  return(orfCodonClassList)
}

#----- Generate ORF Codon Classes -----#

# ~ Read in ORFs from Full Genome ~ #
# Using the same ORF nucleotide positions as for Modified CR6.
fullGenome <- readr::read_file("./CW3_ORFs/CW3_EF014462.1_FullGenome.txt")
fullGenomeNoBreaks <- str_remove_all(fullGenome, pattern = "\\n")

orfSequences <- pmap(.l = list(start = list("ORF1" = 6, 
                                            "ORF2" = 5056,
                                            "ORF3" = 6681,
                                            "ORF4" = 5069),
                               stop = list("ORF1" = 5069,
                                         "ORF2" = 6681,
                                         "ORF3" = 7307,
                                         "ORF4" = 5707)),
                     .f = base::substr,
                     x = fullGenomeNoBreaks)

# save individual ORF sequences
paths <- paste0("./CW3_ORFs/CW3_", names(orfSequences), "_nt.txt")
walk2(.x = orfSequences, .y = paths, 
      .f = ~ (write(x = .x, file = .y)),
      sep = "\t")

# ~ ORF1 ~ #
# Positions: 6 - 5069
# Length: 5064 nt | 1687 aa
orf1CodonClassList <- PopulateCodonClasses(orfFile = "./CW3_ORFs/CW3_ORF1_nt.txt",
                                           orfName = "ORF1",
                                           orfStart = 6, orfEnd = 5069)

orf1CodonLUT <- rep(names(x = orf1CodonClassList), each = 3)
names(orf1CodonLUT) <- seq(from = 6, to = 5069)

orf1 <- new("ORF",
            name = "orf1",
            range = names(orf1CodonLUT),
            LUT = orf1CodonLUT,
            codons = orf1CodonClassList)


# ~ ORF2 ~ #
# Positions 5056 - 6681
# Length: 1626 nt | 541 aa
orf2CodonClassList <- PopulateCodonClasses(orfFile = "./CW3_ORFs/CW3_ORF2_nt.txt",
                                           orfName = "ORF2",
                                           orfStart = 5056, orfEnd = 6681)
orf2CodonLUT <- rep(names(x = orf2CodonClassList), each = 3)
names(orf2CodonLUT) <- seq(from = 5056, to = 6681)

orf2 <- new("ORF",
            name = "orf2",
            range = names(orf2CodonLUT),
            LUT = orf2CodonLUT,
            codons = orf2CodonClassList)

# ~ ORF3 ~ #
# Positions: 6681 - 7307
# Length: 627 nt | 209 aa
orf3CodonClassList <- PopulateCodonClasses(orfFile = "./CW3_ORFs/CW3_ORF3_nt.txt",
                                           orfName = "ORF3",
                                           orfStart = 6681, orfEnd = 7307)

orf3CodonLUT <- rep(names(x = orf3CodonClassList), each = 3)
names(orf3CodonLUT) <- seq(from = 6681, to = 7307)

orf3 <- new("ORF",
            name = "orf3",
            range = names(orf3CodonLUT),
            LUT = orf3CodonLUT,
            codons = orf3CodonClassList)

# ~ ORF4 ~ #
# Positions: 5069 - 5707
# Length: 639 nt | 213 aa

# Make class list and LUT
orf4CodonClassList <- PopulateCodonClasses(orfFile = "./CW3_ORFs/CW3_ORF4_nt.txt",
                                           orfName = "ORF4",
                                           orfStart = 5069, orfEnd = 5707)
orf4CodonLUT <- rep(names(x = orf4CodonClassList), each = 3)
names(orf4CodonLUT) <- seq(from = 5069, to = 5707)

orf4 <- new("ORF",
            name = "orf4",
            range = names(orf4CodonLUT),
            LUT = orf4CodonLUT,
            codons = orf4CodonClassList)

orfList <- list("orf1" = orf1, "orf2" = orf2, "orf3" = orf3, "orf4" = orf4)

# Biostrings AMINO_ACID_CODE with Stop Codon Added
aminoAcidCode <- c(Biostrings::AMINO_ACID_CODE, "*" = "*")

save.image("CW3_ORF_Information.RData")

save(orfList, aminoAcidCode, file = "CW3_ORF_Data.RData")
