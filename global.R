# global.R for variantVieweR_shiny_app

#----- Global variables -----#

previousMtxSize <- 0 # for determining when to re-paint the coverage plot

genome_size <- 7383 # For calculation/formatting avg. genome coverage

#----- Class Definitions -----#

# These classes are populated when the data set is selected

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

#----- Function Definitions -----#

#----- Populating the Sample and Variant Objects -----#

# Returns a list of Sample objects, each of which contains a list of
#   Variant objects.
#dataSetSelect <- "baldridge_rumspringa"
#testDataSet <- GenerateSampleData(dataSetSelect) # used to make dataSetSamples
#dataSetSamples <- as.character(testDataSet$Sample) # send to sampleVector

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
      warning("At least one selected sample has no variants detected. No common variants among all samples.",
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

#----- Set up Sample Data Table -----#

# Alignment Count Data #
GenerateSampleData <- function(dataSet) {
  
  # Read in sample alignment count data, calculate percent MNV
  readCounts <- read.delim(paste0("../", dataSet, "/alignment_counts.txt"))
  readCounts <- readCounts %>%
    mutate("percent_MNV" = round(100*(primary_alignments/total_reads), 
                                 digits = 2))
  
  # Genome Coverage Data #
  # Read in genome coverage count data, calculate avg. genome coverage
  rawFiles <- list.files(path = paste0("../", dataSet, "/sample_data/genome_coverage/"), 
                         pattern = "*_coverage.txt")
  
  # column names determined from documentation at 
  # https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
  columnNames <- c("chromosome", "depth", "number_of_bases",
                   "chromosome_size", "fraction_of_bases")
  
  # Run calculation, store in genomeCovVec - will be combined with readCounts 
  genomeCovVec <- vector(mode = "integer", length = length(rawFiles))
  
  for (i in 1:length(rawFiles)) {
    
    currentFile <- rawFiles[i]
    currentFileName <- as.vector(strsplit(currentFile, "_coverage.txt")[[1]])
    
    currentData <- read.delim(paste0("../", dataSet, 
                                     "/sample_data/genome_coverage/", 
                                     currentFile), 
                              header = FALSE,
                              col.names = columnNames)
    currentData <- currentData %>%
      dplyr::select(depth, number_of_bases) %>%
      dplyr::mutate(product = depth * number_of_bases)
    
    avg_cov <- round(sum(currentData$product)/genome_size,
                     digits = 0)
    
    genomeCovVec[i] <- avg_cov
    names(genomeCovVec)[i] <- currentFileName
    
  }
  
  # Merge Data #
  # Add the average genome coverage numbers to the readCounts data
  sampleData <- readCounts %>%
    dplyr::mutate("avg_genome_cov" = genomeCovVec[sample])
  names(sampleData) <- c("Sample", "Total Reads", "Primary Alignments", 
                         "% MNV", "Average Coverage")
  sampleData$Sample <- as.character(sampleData$Sample)
  
  return(sampleData)
  
}


#----- Generate Sample Variant Table -----#

GetVCF <- function(dataSet, sample) {
  # Read in and format VCF file for selected sample
  variantSample <- paste0("../", dataSet, "/variants/", sample, "_variants.vcf")
  
  vcfFile <- tryCatch({ # in case there are no variants in the VCF file
      read.delim(variantSample,
                 comment.char = "#", # ignore VCF header lines
                 header = FALSE,
                 colClasses = "character") # suppress conversion of columns
    },
    error = function(e) {
      # return an empty data frame
      data.frame("Refence Genome" = character(0), 
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
  
  vcfHeaders <- c("Refence Genome", "Position", 
                  "ID", "Reference", "Alternative",
                  "Quality", "Filter", "Info", "Format", "Values")
  names(vcfFile) <- vcfHeaders
  
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
           "Total Depth" = str_remove(str_extract(Info, "DP=[:digit:]+"),
                                      "DP="),
           # Allelic Freq = (allelic depth/raw depth) * 100%
           "Allelic Frequency" = 
             round((100 * as.numeric(`Allelic Depth`)/as.numeric(`Total Depth`)), 
                                           digits = 2)) %>%
    select(-c("ID", "Filter", "Info", "Format", "Values"))
  
  comment(vcfFileFormatted) <- sample
  
  return(vcfFileFormatted)
  
}

#----- Generate Coverage Plot for Sample & Variants -----#

PlotCoverage <- function(dataSet, sample, positions = NULL, widths = 1) {

  # Get the coverage file
  bedgraphDT <- fread(paste0("../", dataSet, "/alignment_files/", 
                             sample, "_sorted.bedGraph"),
                      col.names = c("chromosome", "start", "end", "value"))

  if (nrow(bedgraphDT) == 0) {
    stop("No coverage data to plot.", call. = FALSE)
  }
  
  # Generate the top axis track
  gtrack <- GenomeAxisTrack(fontsize = 20, fontcolor = "black", col = "black")
  # Generate the coverage track
  dtrack <- DataTrack(range = bedgraphDT, genome = "ModCR6", 
                      type = "histogram", name = " ",
                      background.title = "#2C3E50", col.histogram = "grey28",
                      # hex code matches flatly top bar, previously "slategrey"
                      fontsize = 18)
  
  # is positions null? 
  # if yes - plot tracks w/o highlights 
  # if no - plot tracks with highlights
  if (is.null(positions)) {
    coveragePlot <- plotTracks(list(gtrack, dtrack))
    #trackList <- list("gtrack" = gtrack, "dtrack" = dtrack)
    return(coveragePlot)
  } else {
    htrack <- HighlightTrack(trackList = list(dtrack),
                             start = positions,
                             width= widths,
                             inBackground = FALSE,
                             fill = "#FFE3E6b8")
    coveragePlot <- plotTracks(list(gtrack, htrack))
    #trackList <- list("gtrack" = gtrack, "htrack" = htrack)
    return(coveragePlot)
  }
}

