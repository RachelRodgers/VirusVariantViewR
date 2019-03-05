setClass("Sample",
         representation = representation(sample_name = "character",
                                         variant_list = "list"),
         prototype = prototype(sample_name = NA_character_,
                               variant_list = list()))

setClass("Variant",
         representation = representation(parent_sample = "character",
                                         position = "numeric",
                                         ref_allele = "character",
                                         alt_allele = "character"),
         prototype = prototype(parent_sample = NA_character_,
                               position = NA_real_,
                               ref_allele = NA_character_,
                               alt_allele = NA_character_))

BuildSampleObjects <- function(dataSet, sampleVector) {
  # dataSet will come from input$dataSetSelect
  # sampleVector should automatically populate with all the samples available
  #   in the selected data set
  sampleClassList <- vector(mode = "list", length = length(sampleVector))
  # Build a Sample class object from the vector
  for (i in 1:length(sampleVector)) {
    currentSample <- sampleVector[i]
    # populate a Variant class object
    currentVCFFile <- GetVCF(dataSet, sample = currentSample)
    # Stop here if nothing is in the sample
    if (nrow(currentVCFFile) == 0) {
      warning(paste(currentSample, "has no variants detected"),
              call. = FALSE)
      newSample <- new("Sample",
                       sample_name = currentSample,
                       variant_list = list(0))
      sampleClassList[[i]] <- newSample
      names(sampleClassList)[i] <- currentSample
      #sampleClassList[[currentSample]] <- newSample
    } else {
      currentSampleVariantList <- vector(mode = "list", length = nrow(currentVCFFile))
      for (j in 1:nrow(currentVCFFile)) {
        positionString <- as.character(currentVCFFile[j, "Position"])
        referenceString <- as.character(currentVCFFile[j, "Reference"])
        alternativeString <- currentVCFFile[j, "Alternative"]
        variantID <- paste(positionString, referenceString, alternativeString,
                           sep = "_")
        newVariant <- new("Variant",
                          parent_sample = currentSample,
                          position = as.numeric(positionString),
                          ref_allele = referenceString,
                          alt_allele = alternativeString)
        currentSampleVariantList[[j]] <- newVariant
        names(currentSampleVariantList)[j] <- variantID
        #currentSampleVariantList[[variantID]] <- newVariant
      }
      # Fill in the sample object
      newSample <- new("Sample",
                       sample_name = currentSample,
                       variant_list = currentSampleVariantList)
      # Add to sampleClassList
      sampleClassList[[i]] <- newSample
      names(sampleClassList)[i] <- currentSample
      #sampleClassList[[currentSample]] <- newSample
    }
  } 
  return(sampleClassList)
}

sampleVector <- c("Baldridge_17", "Baldridge_16", "Baldridge_15", "TEST_EMPTY")
dataSet <- "baldridge_rumspringa"

originalSampleObjectList <- BuildSampleObjects(dataSet = dataSet, sampleVector = sampleVector)

FindExactMatches <- function(sampleObjectList) {
  # Check if any samples have 0 variants and remove.  
  variantsAreNull <- vector(mode = "logical", length = length(sampleObjectList))
  for (i in 1:length(sampleObjectList)) {
    variantsAreNull[i]<- is.null(names(sampleObjectList[[i]]@variant_list))
  }
  
  selectedSampleObjectList <- sampleObjectList[!variantsAreNull]
  
  #----- Exact Matching (Position & Variant) -----#
  
  # For all the variants between the selected samples, make a table showing which
  #   samples they are in.
  #   First check that the sample has variants. (is it's variant_list slot == 0?)
  #   If yes, skip it.
  variantDFList <- vector(mode = "list", length = length(selectedSampleObjectList))
  for (k in 1:length(selectedSampleObjectList)) {
    currentSample <- selectedSampleObjectList[[k]]
    currentSampleName <- currentSample@sample_name
    currentSampleDF <- data.frame(names(currentSample@variant_list), 1)
    colnames(currentSampleDF) <- c("variants", currentSampleName)
    variantDFList[[k]] <- currentSampleDF
  }
  
  # Merge the variant data frames together, convert NA's to 0's
  fullDF <- Reduce(f = function(df1, df2) {merge(x = df1, y = df2,
                                                 by = "variants", all = TRUE)},
                   x = variantDFList)
  fullDF[is.na(fullDF)] <- 0
  
  # What are the counts for each variant?
  fullDFModified <- column_to_rownames(fullDF, var = "variants")
  variantCounts <- rowSums(fullDFModified)
  
  # Which samples do you find each variant in?
  variantInSamples <- vector(mode = "character", length = nrow(fullDFModified))
  for (m in 1:nrow(fullDFModified)) {
    currentRow <- fullDFModified[m, ]
    currentVariantName <- rownames(currentRow)
    samplesVecString <- paste(colnames(currentRow)[currentRow != 0], collapse = ", ")
    variantInSamples[m] <- samplesVecString
    names(variantInSamples)[m] <- currentVariantName
  }
  
  # Change things to a DF merge
  variantCountsDF <- data.frame("Number_of_Samples" = variantCounts)
  variantInSamplesDF <- data.frame("Samples" = variantInSamples)
  # Make sure the dimensions of these match
  if (!(identical(dim(variantCountsDF), dim(variantInSamplesDF)))) {
    stop("Something is wrong with the variant/sample counting.",
         call. = FALSE)
  }
  # Otherwise, merge them.
  variantDataDF <- merge(variantCountsDF, variantInSamplesDF, by = "row.names")
  variantDataDF <- variantDataDF %>%
    filter(Number_of_Samples >= 2) %>%
    arrange(desc(Number_of_Samples)) %>%
    column_to_rownames(var = "Row.names")
  
  # Put warning within app.R to deal with empty table.
  
  return(variantDataDF)
}


