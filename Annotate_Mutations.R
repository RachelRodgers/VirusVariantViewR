# Annotate_Mutations.R

# For each VCF file, annotate the reference and alternative mutations.

load("Mod_CR6_ORF_Information.RData")

library("dplyr")
library("stringr")

# Find correct files, store and loop
rawVCFFilePath <- "../larry_mnv_190220/variants"
rawVCFFiles <- list.files(rawVCFFilePath, pattern = ".vcf", full.names = TRUE)
annotatedVCFDirectory <- file.path(rawVCFFilePath, "annotated_variants")
dir.create(annotatedVCFDirectory)

for (j in 1:length(rawVCFFiles)) {
  variantSample <- rawVCFFiles[j]
  variantSampleName <- str_remove(string = basename(variantSample),
                                  pattern = "\\.vcf")
  # for current sample set, open each VCF file, add column for reference and mutant allele
  # if the ref and/or alt is an INDEL, just mark columns as INDEL
  
  vcfFile <- tryCatch({ # in case there are no variants in the VCF file
    read.delim(variantSample,
               comment.char = "#", # ignore VCF header lines
               header = FALSE,
               colClasses = "character") # suppress conversion of columns
  },
  error = function(e) {
    # return an empty data frame
    data.frame("Reference Genome" = character(0), 
               "Position" = character(0), 
               "ID" = character(0), 
               "Reference" = character(0), 
               "Alternative" = character(0),
               "Quality"= character(0), 
               "Filter" = character(0), 
               "Info" = character(0), 
               "Format" = character(0), 
               "Values" = character(0))
  })
  
  vcfHeaders <- c("Reference Genome", "Position", 
                  "ID", "Reference", "Alternative",
                  "Quality", "Filter", "Info", "Format", "Values")
  names(vcfFile) <- vcfHeaders
  
  # If VCF file isn't empty, add a column for Reference Protein, Alternative Protein,
  #   and Mutation Type
  
  # set orf ranges
  orf1Range <- seq(from = 6, to = 5069)
  orf2Range <- seq(from = 5070, to = 6681)
  orf3Range <- seq(from = 6681, to = 7307)
  
  if (nrow(vcfFile) != 0) {
    vcfFileAnnotated <- cbind(vcfFile, 
                              "Reference Codon" = NA, "Reference Protein" = NA,
                              "Mutant Codon" = NA, "Mutant Protein" = NA, 
                              "Mutation Type" = NA)
    
    # loop over the rows and populate the new columns
    for (i in 1:nrow(vcfFileAnnotated)) {
      
      currentRow <- vcfFileAnnotated[i, ] 
      
      if (nchar(currentRow$Reference) > 1 | nchar(currentRow$Alternative) > 1) {
        vcfFileAnnotated[i, "Mutation Type"] <- "INDEL"
        next()
        
      } else {
        mutantAllele <- currentRow$Alternative
        mutantPosition <- currentRow$Position
        
        currentCodonLUT <- NULL
        currentCodonClassList <- NULL
        
        if (mutantPosition %in% orf1Range) {
          currentCodonLUT <- orf1CodonLUT
          currentCodonClassList <- orf1CodonClassList
        } else if (mutantPosition %in% orf2Range) {
          currentCodonLUT <- orf2CodonLUT
          currentCodonClassList <- orf2CodonClassList
        } else if (mutantPosition %in% orf3Range) {
          currentCodonLUT <- orf3CodonLUT
          currentCodonClassList <- orf3CodonClassList
        }
        
        # if mutant position doesn't align to any orf, move on to next mutation
        if (is.null(currentCodonClassList)) {
          currentRow$`Mutation Type` <- "non-coding"
          isNonCodingMutation <- TRUE
          next()
        }
        
        # determine identity of ref protein
        referenceCodonName <- currentCodonLUT[as.character(mutantPosition)]
        referenceCodon <- currentCodonClassList[[referenceCodonName]]
        referenceCodonString <- paste(referenceCodon@sequence_vector, collapse = "")
        referenceProtein <- Biostrings::GENETIC_CODE[[referenceCodonString]]
        referenceProteinFull <- aminoAcidCode[[referenceProtein]]
        vcfFileAnnotated[i, "Reference Protein"] <- referenceProteinFull
        vcfFileAnnotated[i, "Reference Codon"] <- referenceCodonString
        
        # determine identity of mutant protein from ref codon
        mutantCodon <- referenceCodon@sequence_vector
        mutantCodon[as.character(mutantPosition)] <- mutantAllele
        mutantCodonString <- paste(mutantCodon, collapse = "")
        mutantProtein <- Biostrings::GENETIC_CODE[[mutantCodonString]]
        mutantProteinFull <- aminoAcidCode[[mutantProtein]]
        vcfFileAnnotated[i, "Mutant Protein"] <- mutantProteinFull
        vcfFileAnnotated[i, "Mutant Codon"] <- mutantCodonString
        
        # Determine mutation type
        if (referenceProtein == mutantProtein) {
          # synonymous
          vcfFileAnnotated[i, "Mutation Type"] <- "synonymous"
        } else {
          if (mutantProtein == "*") {
            # nonsense
            vcfFileAnnotated[i, "Mutation Type"] <- "nonsense"
          } else {
            # missense
            vcfFileAnnotated[i, "Mutation Type"] <- "missense"
          }
        }
        
      }
    }
    
    # write the file out
    write.table(x = vcfFileAnnotated,
                file = paste(annotatedVCFDirectory, "/", 
                             variantSampleName, "_annotated.txt",
                             sep = ""),
                append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
    
  } else { # empty file
    write.table(x = vcfFile,
                file = paste(annotatedVCFDirectory, "/",
                             variantSampleName, "_annotated.txt",
                             sep = ""),
                append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  
  
}

