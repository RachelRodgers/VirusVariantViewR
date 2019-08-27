# GetVCF.R

source("./global.R")

#----- Generate Sample Variant Table -----#
dataSet <- "Combined_Data"
sampleData <- GenerateSampleData(dataSet = dataSet)
sampleList <- as.character(sampleData$Sample)

# Find correct files, store and loop
annotatedVCFFilePath <- paste0("../", dataSet, "/variants/annotated_variants")
annotatedVCFFiles <- list.files(rawVCFFilePath, pattern = ".vcf", full.names = TRUE)
formattedVCFDirectory <- file.path(annotatedVCFFilePath, "formatted_variants")
dir.create(formattedVCFDirectory)

for (i in 1:length(sampleList))
  sample <- sampleList[1]
  # Read in and format VCF file for selected sample
  vcfFile <- data.table::fread(file = paste0("../", dataSet, 
                                             "/variants/annotated_variants/", 
                                             sample, "_variants_annotated.txt"),
                               header = TRUE, sep = "\t", 
                               stringsAsFactors = FALSE,
                               verbose = FALSE)
  
  # Read in the bedGraphDT file to get Total Depth Information, append to vcfFile
  bedgraphDT <- fread(paste0("../", dataSet, "/alignment_files/",
                             sample, "_sorted.bedGraph"),
                      col.names = c("chromosome", "start", "end", "value"),
                      colClasses = c("character", "integer", "integer", "character"))
  
  coverageRangeList <- vector(mode = "list", length = nrow(bedgraphDT))
  
  for (i in 1:nrow(bedgraphDT)) {
    coverageRangeList[[i]] <- seq(from = as.numeric(bedgraphDT[i, "start"]),
                                  to = as.numeric(bedgraphDT[i, "end"] - 1)) # -1 because end is also start of next range
    names(coverageRangeList)[i] <- bedgraphDT[i, "value"]
  }
  
  GetCoverageAtPosition <- function(position) {
    for (i in 1:length(coverageRangeList)) {
      if (position %in% coverageRangeList[[i]]) {
        covValue <- names(coverageRangeList)[i]
        return(covValue)
        next()
      }
    }
  }

  for (i in 1:nrow(vcfFile)) {
    currentPosition <- vcfFile[i, "Position"]
    vcfFile[i, "Total Depth"] <- GetCoverageAtPosition(position = currentPosition)
  }
  
  # Get the primary alignment value from the current sample's alignment counts
  sampleAlignmentCounts <- GenerateSampleData(dataSet)
  samplePrimaryAlignments <- subset(sampleAlignmentCounts,
                                    Sample == sample)$`Primary Alignments`
  # Filter sampleAlignmentCounts by sample
  #samplePrimaryAlignments <- subset(GenerateSampleData,
  #sample == sample)$primary_alignments
  
  # The raw read depth at each position will be given by DP=xx; the first element
  #   in the INFO field.  Will divide the allelic depth by this value to get
  #   the allelic frequency.
  # The last number in Values will be the allelic depth - unfiltered number of 
  #   reads supporting the reported allele(s)
  
  vcfFileFormatted <- vcfFile %>%
    mutate("Allelic Depth" = map_chr(.x = str_split(Values, pattern = ","),
                                     .f = tail, n = 1),
           #"Total Depth" = str_remove(str_extract(Info, "DP=[:digit:]+"), "DP="),
           # Allelic Freq = (allelic depth/raw depth) * 100%
           "Allelic Frequency (%)" = 
             round((100 * as.numeric(`Allelic Depth`)/as.numeric(`Total Depth`)), 
                   digits = 2)) %>%
    select(-c("ID", "Filter", "Info", "Format", "Values"))
  
  #write.file()
  
  #comment(vcfFileFormatted) <- sample
  
  #return(vcfFileFormatted)
  
  
  
}