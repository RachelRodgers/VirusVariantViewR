# StaticCoveragePlots.R

# Script used to generate static coverage plots of the Combined_Data samples.
#   For each sample, generate the gtrack, the dtrack (requires bedGraph file),
#   and the htrack (requires VCF file).

source("./global.R")
options(ucscChromosomeNames = FALSE)

setClass("UpdatedSample",
         representation = representation(no_variants = "logical"),
         contains = "Sample")

setClass("VariantCoverage",
         representation = representation(DP = "character"),
         contains = "Variant")

# Returns a list of Sample objects, each of which contains a list of
#   Variant objects.
dataSet <- "Combined_Data"
sampleData <- GenerateSampleData(dataSet) # used to make dataSetSamples
sampleVector <- as.character(sampleData$Sample) # send to sampleVector

BuildSampleObjectsModified <- function(dataSet, sampleVector) {
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
      newSample <- new("UpdatedSample",
                       sample_name = currentSample,
                       variant_list = list(0),
                       no_variants = TRUE)
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
        totalDepth <- as.character(currentVCFFile[j, "Total Depth"])
        variantDF <- data.frame("Variation" = variantID,
                                "Variation_Readable" = variantReadableID,
                                "Reference_Allele" = referenceString,
                                "Reference_Codon" = currentVCFFile[j, "Reference Codon"],
                                "Reference_Protein" = currentVCFFile[j, "Reference Protein"],
                                "Mutant_Allele" = alternativeString,
                                "Mutant_Codon" = currentVCFFile[j, "Mutant Codon"],
                                "Mutant_Protein" = currentVCFFile[j, "Mutant Protein"],
                                "Mutation_Type" = currentVCFFile[j, "Mutation Type"])
        newVariant <- new("VariantCoverage",
                          parent_sample = currentSample,
                          position = as.numeric(positionString),
                          ref_allele = referenceString,
                          alt_allele = alternativeString,
                          variant_df = variantDF,
                          DP = totalDepth)
        currentSampleVariantList[[j]] <- newVariant
        names(currentSampleVariantList)[j] <- variantID
      }
      
      # Fill in the sample object
      newSample <- new("UpdatedSample",
                       sample_name = currentSample,
                       variant_list = currentSampleVariantList,
                       no_variants = FALSE)
      # Add to sampleClassList
      sampleClassList[[i]] <- newSample
      names(sampleClassList)[i] <- currentSample
      #sampleClassList[[currentSample]] <- newSample
    }
  } 
  return(sampleClassList)
}

sampleObjectsList <- BuildSampleObjectsModified(dataSet = dataSet,
                                                sampleVector = sampleVector)
for (i in 1:length(sampleObjectsList)) {
  currentSample <- sampleObjectsList[[i]]
  currentSampleName <- currentSample@sample_name
  
  bedgraphDT <- fread(paste0("../", dataSet, "/alignment_files/", 
                             currentSampleName, "_sorted.bedGraph"),
                      col.names = c("chromosome", "start", "end", "value"))
  # Check if there's any coverage data to plot, if not, go to next
  if (nrow(bedgraphDT) == 0) {
    stop("No coverage data to plot.", call. = FALSE)
    next()
  }
  
  # Otherwise, there is coverage data, so get file path ready
  plotPath <- paste0("../Combined_Data/static_coverage_plots/",
                     currentSampleName, ".pdf")
  
  # Now check if there are any variants to loop over
  noVariants <- currentSample@no_variants
  
  # If there are variants to look at, loop over them.  allPositionsNum will
  #   hold the passing variants, if there are any.  Otherwise it'll be NULL.
  if (noVariants == FALSE) {
    currentVariantList <- currentSample@variant_list
    allPositionsNum <- NULL
    allPositions <- vector(mode = "character")
    
    for (j in 1:length(currentVariantList)) {
      # first check that depth is greater than 100
      currentTotalDepth <- as.numeric(currentVariantList[[j]]@DP)
      if(currentTotalDepth >= 100) {
        allPositions <- c(allPositions, currentVariantList[[j]]@position)
      }
    }
    
    allPositionsNum <- as.numeric(allPositions)
  }
  
  # Generate the genome  and data track, which will exist whether there are 
  #   variants to highlight or not
  gtrack <- GenomeAxisTrack(fontsize = 20, fontcolor = "black", col = "black")
  
  dtrack <- DataTrack(range = bedgraphDT, genome = "ModCR6", 
                      type = "histogram", name = " ",
                      background.title = "#2C3E50", col.histogram = "grey47",
                      fontsize = 18,
                      baseline = 100, lty.baseline = 2, lwd.baseline = 2,
                      col.baseline = "black")
  
  if (noVariants == TRUE | is.null(allPositionsNum) | length(allPositionsNum) == 0) {
    pdf(file = plotPath, paper = "USr", width = 9, height = 4)
    coveragePlot <- plotTracks(list(gtrack, dtrack))
    dev.off()
  } else {
    # Generate the highlight track
    htrack <- HighlightTrack(trackList = list(dtrack),
                             start = allPositionsNum,
                             width = 1,
                             inBackground = FALSE,
                             fill = "#FFE3E6b8")
    
    pdf(file = plotPath, paper = "USr", width = 9, height = 4)
    coveragePlot <- plotTracks(list(gtrack, htrack))
    dev.off()
  }
}
