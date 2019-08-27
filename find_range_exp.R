source("./global.R")

dataSet<- "Combined_Data"
testDataSet <- GenerateSampleData(dataSet) # used to make dataSetSamples
dataSetSamples <- as.character(testDataSet$Sample) # send to sampleVector

sampleObjectList <- BuildSampleObjects(dataSet = dataSet,
                                       sampleVector = dataSetSamples)

sample <- sampleObjectList[[1]]@sample_name

# w/n GetVCF()
vcfFile <- data.table::fread(file = paste0("../", dataSet, 
                                           "/variants/annotated_variants/", 
                                           sample, "_variants_annotated.txt"),
                             header = TRUE, sep = "\t", 
                             stringsAsFactors = FALSE,
                             verbose = FALSE)

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

for (j in 1:length(vcfFile)) {
  currentPosition <- vcfFile[j, "Position"]
  vcfFile[j, "Total Depth"] <- GetCoverageAtPosition(position = currentPosition)
}


