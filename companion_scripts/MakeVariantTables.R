# MakeVariantTables.R

# This script is used to make a presence/absence table of variants for every
#   sample in the specified data set. It should also contain any metadata available for
#   the samples.  

# I will be starting with the rumspringa "Combined_Data" set of samples, where
#   there's no explicit metadata.

source("global.R")

# Grab the sample data
dataSetName <- "Combined_Data"
sampleData <- GenerateSampleData(dataSet = dataSetName)
dataSetSamples <- as.character(sampleData$Sample)

# Build the sample class objects
sampleObjectList <- BuildSampleObjects(dataSet = dataSetName,
                                       sampleVector = dataSetSamples)

# Combine each samples' variants into one large data frame:

# First check if any samples have 0 variants and remove from sampleObjectList.  
variantsAreNull <- vector(mode = "logical", length = length(sampleObjectList))
for (i in 1:length(sampleObjectList)) {
  variantsAreNull[i]<- is.null(names(sampleObjectList[[i]]@variant_list))
}

whosNull <- names(sampleObjectList[variantsAreNull])
whosNull

selectedSampleObjectList <- sampleObjectList[!variantsAreNull]

# Stop if all the samples selected have no variants in them.
if (length(selectedSampleObjectList) == 0) {
  stop("No variants are detected in the selected samples.  Cannot identify common variants.",
       call. = FALSE)
}

# For all the variants between the selected samples, make a table showing which
#   samples they are in.

sampleVariantsList <- vector(mode = "list", length = length(selectedSampleObjectList))

for (k in 1:length(selectedSampleObjectList)) {
  currentSample <- selectedSampleObjectList[[k]]
  currentSampleDF <- data.frame(names(currentSample@variant_list), 1)
  colnames(currentSampleDF) <- c("variants", currentSample@sample_name)
  sampleVariantsList[[k]] <- currentSampleDF
}

combinedVariantDF <- Reduce(f = function(df1, df2) {merge(x = df1, y = df2,
                                                          by = "variants",
                                                          all = TRUE)},
                             x = sampleVariantsList)

combinedVariantDF[is.na(combinedVariantDF)] <- 0

combinedVariantDF$variants <- as.character(combinedVariantDF$variants)
arrVariantDF <- arrange(combinedVariantDF, variants)

write.table(arrVariantDF, file = "../Combined_Data/VariantPATable.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
