# Annotate_Mutations.R

# For each VCF file, annotate the reference and alternative mutations and add
#   total depth information from the sample's bedGraph file.

load("ORF_Data.RData")

library("dplyr")
library("stringr")
library("data.table")

# Find correct files, store and loop
dataSet <- "baldridge_rumspringa"
rawVCFFilePath <- paste0("../", dataSet, "/variants")
rawVCFFiles <- list.files(rawVCFFilePath, pattern = ".vcf", full.names = TRUE)
annotatedVCFDirectory <- file.path(rawVCFFilePath, "annotated_variants")
dir.create(annotatedVCFDirectory)

for (i in 1:length(rawVCFFiles)) {
  variantSample <- rawVCFFiles[i]
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
  if (nrow(vcfFile) == 0) {
    write.table(x = vcfFile,
                file = paste(annotatedVCFDirectory, "/",
                             variantSampleName, "_annotated.txt",
                             sep = ""),
                append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  } else {
    vcfFileAnnotated <- cbind(vcfFile, "Location" = NA,
                              "Reference Codon" = NA, "Reference Protein" = NA,
                              "Mutant Codon" = NA, "Mutant Protein" = NA, 
                              "Mutation Type" = NA)
    
    # loop over the rows and populate the new columns
    for (j in 1:nrow(vcfFileAnnotated)) {
      mutationIsINDEL <- FALSE
      currentRow <- vcfFileAnnotated[j, ] 
      # Mutation Information:
      mutantPosition <- currentRow$Position
      mutantAllele <- currentRow$Alternative
      
      # Do we have an INDEL?
      if (nchar(currentRow$Reference) > 1 | nchar(currentRow$Alternative) > 1) {
        mutationIsINDEL <- TRUE
      } 
      
      # Which orf(s) is this position in?
      orfSelectionVec <- vector(mode = "logical", length = length(orfList))
      for (k in 1:length(orfList)) {
        orfSelectionVec[k] <- mutantPosition %in% orfList[[k]]@range
      }
      
      # Get the correct orf objects
      currentORFObjects <- orfList[orfSelectionVec]
      orfCount <- length(currentORFObjects)
      
      # In case no ORF objects were selected:
      if (orfCount == 0) {
        vcfFileAnnotated[j, "Location"] <- "Non-coding"
        next() # go to next line in file
      }
      
      # Loop over the ORF objects and extract the location information
      # Make vector to hold location information in case there is more than one orf.
      #   This will be concatenated together later and appear in a single cell in the
      #   variant table.
      locationInfoVec <- vector("character", orfCount)
      referenceCodonVec <- vector("character", orfCount)
      referenceProteinVec <- vector("character", orfCount)
      mutantCodonVec <- vector("character", orfCount)
      mutantProteinVec <- vector("character", orfCount)
      mutationTypeVec <- vector("character", orfCount)
      
      for (k in 1:length(currentORFObjects)) {
        # What's the current ORF?
        currentORF <- currentORFObjects[[k]]
        currentORFName <- currentORF@name
        currentORFCodons <- currentORF@codons
        # What codon should we focus on?
        referenceCodonName <- currentORF@LUT[mutantPosition]
        # ~ Location Info ~ #
        currentGene <- currentORFCodons[[referenceCodonName]]@gene_name
        aaPosition <- currentORFCodons[[referenceCodonName]]@protein_position
        locationString <- paste(toupper(currentORFName), currentGene, "AA#", aaPosition)
        locationInfoVec[k] <- locationString
        
        # ~ Is mutation an INDEL? ~ #
        if (mutationIsINDEL == TRUE) {
          vcfFileAnnotated[j, "Mutation Type"] <- "INDEL"
          break() # exit this loop
        }
        
        # ~ Reference Information ~ #
        referenceCodonString <- paste(currentORFCodons[[referenceCodonName]]@sequence_vector, 
                                      collapse = "")
        referenceCodonVec[k] <- referenceCodonString
        
        referenceProtein <- aminoAcidCode[[Biostrings::GENETIC_CODE[[referenceCodonString]]]]
        referenceProteinVec[k] <- referenceProtein
        
        # ~ Mutant Information ~ #
        mutantCodon <- currentORFCodons[[referenceCodonName]]@sequence_vector
        mutantCodon[as.character(mutantPosition)] <- mutantAllele # stick in the mutation
        mutantCodonString <- paste(mutantCodon, collapse = "")
        mutantCodonVec[k] <- mutantCodonString
        
        mutantProtein <- aminoAcidCode[[Biostrings::GENETIC_CODE[[mutantCodonString]]]]
        mutantProteinVec[k] <- mutantProtein
        
        # ~ Mutation Type ~ #
        mutationTypeVec[k] <- case_when(referenceProtein == mutantProtein ~ "synonymous",
                                        mutantProtein == "*" ~ "nonsense",
                                        TRUE ~ "missense")
      } # done looping over ORF objects 
      
      # Populate the annotations.
      vcfFileAnnotated[j, "Location"] <- paste(locationInfoVec, collapse = "; ")
      
      if (mutationIsINDEL == TRUE) {
        next() # go to next variant
      }
      
      vcfFileAnnotated[j, "Reference Codon"] <- paste(referenceCodonVec, collapse = "; ")
      vcfFileAnnotated[j, "Reference Protein"] <- paste(referenceProteinVec, collapse = "; ")
      vcfFileAnnotated[j, "Mutant Codon"] <- paste(mutantCodonVec, collapse = "; ")
      vcfFileAnnotated[j, "Mutant Protein"] <- paste(mutantProteinVec, collapse = "; ")
      vcfFileAnnotated[j, "Mutation Type"] <- paste(mutationTypeVec, collapse = "; ")
    }
    
    # Add Total Depth information ere, then write out (see find_range_exp.R)
    # Get sample Name
    sampleName <- str_remove(variantSampleName, "_variants")
    # read in bedgraphDT
    bedgraphDT <- fread(paste0("../", dataSet, "/alignment_files/",
                               sampleName, "_sorted.bedGraph"),
                        col.names = c("chromosome", "start", "end", "value"),
                        colClasses = c("character", "integer", "integer", "character"))
    # generate coverage map
    coverageRangeList <- vector(mode = "list", length = nrow(bedgraphDT))
    
    for (j in 1:nrow(bedgraphDT)) {
      coverageRangeList[[j]] <- seq(from = as.numeric(bedgraphDT[j, "start"]),
                                    to = as.numeric(bedgraphDT[j, "end"] - 1)) # -1 because end is also start of next range
      names(coverageRangeList)[j] <- bedgraphDT[j, "value"]
    }
    
    GetCoverageAtPosition <- function(position) {
      for (pos_idx in 1:length(coverageRangeList)) {
        if (position %in% coverageRangeList[[pos_idx]]) {
          covValue <- names(coverageRangeList)[pos_idx]
          return(covValue)
          next()
        }
      }
    }
    
    for (j in 1:nrow(vcfFileAnnotated)) {
      currentPosition <- vcfFileAnnotated[j, "Position"]
      vcfFileAnnotated[j, "Total Depth"] <- GetCoverageAtPosition(position = currentPosition)
    }
    
    # Calculate Allelic Depth
    
    # The AD value is given as ##,## and is the last key-value element in the
    #   Value column of the vcf file.  The leading ## represents the number of
    #   reads supporting the reference allele, while the second ## represents
    #   the number of reads supporting the alternative allele.  Dividing
    #   the alt reads by (alt + ref reads) should yield the allelic frequency.
    
    vcfFileAnnotated <- vcfFileAnnotated %>% 
      mutate("AD" = map_chr(.x = str_split(Values, pattern = ":"), .f = tail, n = 1),
             "Ref_Reads" = map_chr(.x = str_split(AD, pattern = ","), .f = head, n = 1),
             "Alt_Reads" = map_chr(.x = str_split(AD, pattern = ","), .f = tail, n = 1),
             "Allelic Frequencey (%)" = ifelse(Ref_Reads == "0",
                                               yes = "100",
                                               no = 
                                                 round(100 * (as.numeric(Alt_Reads)/(as.numeric(Ref_Reads) + as.numeric(Alt_Reads))),
                                                       digits = 2)))
    
    # write the file out
    write.table(x = vcfFileAnnotated,
                file = paste(annotatedVCFDirectory, "/",
                             variantSampleName, "_annotated.txt",
                             sep = ""),
                append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  }
}
