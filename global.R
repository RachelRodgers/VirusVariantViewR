# global.R for variantVieweR_shiny_app

#----- Global variables -----#

previousMtxSize <- 0 # for determining when to re-paint the coverage plot

genome_size <- 7383 # For calculation/formatting avg. genome coverage

#----- Function Definitions -----#

#----- Set up Sample Data Table -----#

# Alignment Count Data #
GenerateAlignmentCounts <- function(dataSet) {
  
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
  
  vcfFile <- read.delim(variantSample,
                        comment.char = "#", # ignore VCF header lines
                        header = FALSE,
                        colClasses = "character") # suppress conversion of columns
  vcfHeaders <- c("Refence Genome", "Position", 
                  "ID", "Reference", "Alternative",
                  "Quality", "Filter", "Info", "Format", "Values")
  names(vcfFile) <- vcfHeaders
  
  # Get the primary alignment value from the current sample's alignment counts
  sampleAlignmentCounts <- GenerateAlignmentCounts(dataSet)
  samplePrimaryAlignments <- subset(sampleAlignmentCounts,
                                    Sample == sample)$`Primary Alignments`
  # Filter sampleAlignmentCounts by sample
  #samplePrimaryAlignments <- subset(GenerateAlignmentCounts,
                                    #sample == sample)$primary_alignments
  
  # The last number in Values will be the allelic depth - unfiltered number of 
  #   reads supporting the reported allele(s)
  vcfFileFormatted <- vcfFile %>%
    mutate("Allelic Depth" = map_chr(.x = str_split(Values, pattern = ","),
                                     .f = tail, n = 1),
           "Allelic Frequency" = round((100 * 
                                         as.numeric(`Allelic Depth`)/samplePrimaryAlignments),
                                       digits = 2)) %>%
    select(-c("ID", "Filter", "Info", "Format", "Values"))
}

#----- Generate Coverage Plot for Sample & Variants -----#

PlotCoverage <- function(dataSet, sample, positions = NULL, widths = 1) {

  # Get the coverage file
  bedgraphDT <- fread(paste0("../", dataSet, "/alignment_files/",
                             sample, "_sorted.bedGraph"),
                      col.names = c("chromosome", "start", "end", "value"))
  # Generate the top axis track
  gtrack <- GenomeAxisTrack(fontsize = 20, fontcolor = "black", col = "black")
  # Generate the coverage track
  dtrack <- DataTrack(range = bedgraphDT, genome = "ModCR6", 
                      type = "histogram", name = "Coverage",
                      background.title = "#2C3E50", col.histogram = "grey28",
                      # hex code matches flatly top bar, previously "slategrey"
                      fontsize = 20)
  
  # is positions null? 
  # if yes - plot tracks w/o highlights 
  # if no - plot tracks with highlights
  if (is.null(positions)) {
    coveragePlot <- plotTracks(list(gtrack, dtrack))
  } else {
    htrack <- HighlightTrack(trackList = list(dtrack),
                             start = positions,
                             width= widths,
                             inBackground = FALSE,
                             fill = "#FFE3E6b8")
    coveragePlot <- plotTracks(list(gtrack, htrack))
  }
}

