setClass("Sample",
         representation = representation(sample_name = "character",
                                         variant_list = "list"),
         prototype = prototype(sample_name = NA_character_,
                               variant_list = list()))

setClass("Variant",
         representation = representation(parent_sample = "character",
                                         position = "numeric",
                                         ref_allele = "character",
                                         alt_allele = "character",
                                         variant_df = "data.frame"))

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
    } else {
      currentSampleVariantList <- vector(mode = "list", length = nrow(currentVCFFile))
      
      for (j in 1:nrow(currentVCFFile)) {
        positionString <- as.character(currentVCFFile[j, "Position"])
        referenceString <- as.character(currentVCFFile[j, "Reference"])
        alternativeString <- currentVCFFile[j, "Alternative"]
        variantID <- paste(positionString, referenceString, alternativeString,
                           sep = "_")
        variantReadableID <- paste(positionString, " ", referenceString, "/", alternativeString, sep = "")
        variantDF <- data.frame("Variation" = variantID,
                                "Variation_Readable" = variantReadableID,
                                "Reference_Allele" = referenceString,
                                "Reference_Codon" = currentVCFFile[j, "Reference Codon"],
                                "Reference_Protein" = currentVCFFile[j, "Reference Protein"],
                                "Mutant_Allele" = alternativeString,
                                "Mutant_Codon" = currentVCFFile[j, "Mutant Codon"],
                                "Mutant_Protein" = currentVCFFile[j, "Mutant Protein"],
                                "Mutation_Type" = currentVCFFile[j, "Mutation Type"])
        newVariant <- new("Variant",
                          parent_sample = currentSample,
                          position = as.numeric(positionString),
                          ref_allele = referenceString,
                          alt_allele = alternativeString,
                          variant_df = variantDF)
        currentSampleVariantList[[j]] <- newVariant
        names(currentSampleVariantList)[j] <- variantID
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


sampleObjectList <- originalSampleObjectList
FindExactMatches <- function(sampleObjectList) {
  # Check if any samples have 0 variants and remove.  
  variantsAreNull <- vector(mode = "logical", length = length(sampleObjectList))
  for (i in 1:length(sampleObjectList)) {
    variantsAreNull[i]<- is.null(names(sampleObjectList[[i]]@variant_list))
  }
  
  selectedSampleObjectList <- sampleObjectList[!variantsAreNull]
  
  # Stop if all the samples selected have no variants in them.
  if (length(selectedSampleObjectList) == 0) {
    stop("No variants are detected in the selected samples.  Cannot identify common variants.",
         call. = FALSE)
  }
  
  #----- Exact Matching (Position & Variant) -----#
  
  # For all the variants between the selected samples, make a table showing which
  #   samples they are in.
  
  sampleVariantsDFList <- vector(mode = "list", length = length(selectedSampleObjectList))
  sampleVariantsPresenceList <- vector(mode = "list", length = length(selectedSampleObjectList))
  # pull extra info here and stick into currentSampleDF
  for (k in 1:length(selectedSampleObjectList)) {
    # Get current sample information
    currentSample <- selectedSampleObjectList[[k]]
    currentSampleName <- currentSample@sample_name
    
    # Build concatenated data frame from all the variants' variant_df slots
    #   get current Variants, put in list
    currentVariantList <- currentSample@variant_list
    #   pull out each variant's variant_df, put in list
    currentVariantDFList <- map(.x = seq(1:length(currentVariantList)),
                                .f = function(idx) currentVariantList[[idx]]@variant_df)
    #   reduce the list into one big df - (all variants' info for one sample only)
    combinedVariantDFs <- Reduce(f = function(df1, df2) {rbind(x = df1, y = df2)},
                                 x = currentVariantDFList)
    
    #   store the samples variants info in sampleVariantsDFList to be reduced later
    sampleVariantsDFList[[k]] <- combinedVariantDFs
    
    # Build variantPresenceDF that will be used for searching for common variants
    variantPresenceDF <- data.frame(names(currentSample@variant_list), 1)
    colnames(variantPresenceDF) <- c("variants", currentSampleName)
    sampleVariantsPresenceList[[k]] <- variantPresenceDF
  }
  
  # Merge the presence data frames together, convert NA's to 0's
  fullPresenceDF <- Reduce(f = function(df1, df2) {merge(x = df1, y = df2,
                                                 by = "variants", all = TRUE)},
                           x = sampleVariantsPresenceList)
  fullPresenceDF[is.na(fullPresenceDF)] <- 0
  
  # remove extra columns before attempting rowSums
  
  # What are the counts for each variant?
  fullPresenceDFModified <- column_to_rownames(fullPresenceDF, var = "variants")
  variantCounts <- rowSums(fullPresenceDFModified)
  
  # Which samples do you find each variant in?
  variantInSamples <- vector(mode = "character", length = nrow(fullPresenceDFModified))
  for (m in 1:nrow(fullPresenceDFModified)) {
    currentRow <- fullPresenceDFModified[m, ]
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
    column_to_rownames(var = "Row.names") %>%
    rownames_to_column("Variation")
  
  # Put warning within app.R to deal with empty table.
  
  # Now pull extra information before displaying table - from Reducing information 
  #   from sampleVariantsDFList
  variantInformationDF <- Reduce(f = function(df1, df2) {merge(x = df1, y = df2, all = TRUE)},
                                 x = sampleVariantsDFList)
  
  # From variantInformationDF, append information to variantDataDF
  variantInformationFullDF <- merge(x = variantDataDF, variantInformationDF, 
                                    by = "Variation", all = FALSE)
  
  # Modify the table
  variantInformationFinal <- variantInformationFullDF %>%
    select(-c(Variation)) %>%
    select(Variation_Readable, everything()) %>%
    data.table::setnames(old = "Variation_Readable",
                         new = "Variation")
  
  return(variantInformationFinal)
}

## TEST ##
sampleVector <- c("Baldridge_Filtered", "Baldridge_Unfiltered")
dataSet <- "larry_mnv_190206"

originalSampleObjectList <- BuildSampleObjects(dataSet = dataSet, sampleVector = sampleVector)
matchesTest <- FindExactMatches(sampleObjectList = originalSampleObjectList)
