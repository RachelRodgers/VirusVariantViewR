#----- Libraries -----#

library("shiny")
library("DT")
library("Gviz")
library("data.table")
library("tidyverse")

options(shiny.sanitize.errors = FALSE) # need to see the error
options(ucscChromosomeNames = FALSE) # for Gvis

#----- Global variables -----#

previousMtxSize <- 0 # for determining when to re-paint the coverage plot

genome_size <- 7383 # For calculation/formatting avg. genome coverage

#----- Function Definitions -----#

GetVCF <- function(sample) {
  # Read in and format VCF file for selected sample
  variantSample <- paste0("../variants/", sample, "_variants.vcf")

  vcfFile <- read.delim(variantSample,
                        comment.char = "#", # ignore VCF header lines
                        header = FALSE,
                        colClasses = "character") # suppress conversion of columns
  vcfHeaders <- c("Refence Genome", "Position", 
                  "ID", "Reference", "Alternative",
                  "Quality", "Filter", "Info", "Format", "Values")
  names(vcfFile) <- vcfHeaders
  # The last number in Values will be the allelic depth - unfiltered number of 
  #   reads supporting the reported allele(s)
  vcfFileFormatted <- vcfFile %>%
    mutate("Allelic Depth" = map_chr(.x = str_split(Values, pattern = ","),
                                    .f = tail, n = 1)) %>%
    select(-c("ID", "Filter", "Info", "Format", "Values"))
}

PlotCoverage <- function(sample, positions = NULL, widths = 1) {
  
  # Get the coverage file
  bedgraphDT <- fread(paste0("../alignment_files/", sample, "_sorted.bedGraph"),
                      col.names = c("chromosome", "start", "end", "value"))
  # Generate the top axis track
  gtrack <- GenomeAxisTrack(fontsize = 20, fontcolor = "black", col = "black")
  # Generate the coverage track
  dtrack <- DataTrack(range = bedgraphDT, genome = "ModCR6", 
                      type = "histogram", name = "Coverage",
                      background.title = "slategrey", col.histogram = "grey28",
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

#----- Set up Data Table for App -----#

# Alignment Count Data #
# Read in sample alignment count data, calculate percent MNV
readCounts <- read.delim("../alignment_counts.txt")
readCounts <- readCounts %>%
  mutate("percent_MNV" = round(100*(primary_alignments/total_reads), 
                               digits = 2))

# Genome Coverage Data #
# Read in genome coverage count data, calculate avg. genome coverage
rawFiles <- list.files(path = "../sample_data/genome_coverage/", 
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
  
  currentData <- read.delim(paste0("../sample_data/genome_coverage/", 
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


#--------------------- Format, Run App ---------------------#

ui <- navbarPage(
  
  # title id refers to navbarPage object
  title = "MNV (CR6) Variants", id = "navbarPage",
  
  tabPanel(title = "Table", DT::dataTableOutput("sampleDataTable")),
  
  tabPanel(title = "Variants", 
           textInput("varTabInput", "Row ID"),
           DT::DTOutput("varTable"),
           plotOutput("coveragePlot"))
)

server <- function(input, output, session) {
  
  # columnDefs - hide the 3rd column (number primary alignments) 
  #   which is index 2 in the DT object
  # formatStyle - add CSS style 'cursor: pointer' to the 1st column (i.e. sample)
  #   which is index 1 in in the DT object
  output$sampleDataTable <- DT::renderDataTable({
    
    data = datatable(sampleData, selection = "single", rownames = FALSE,
                     options = list(
                       columnDefs = list(list(visible = FALSE, 
                                              targets = c(2))))) %>%
      formatStyle(columns = 1, cursor = "pointer") 
    
  })
  
  # see https://yihui.shinyapps.io/DT-click/
  
  observeEvent(input$sampleDataTable_cell_clicked, {
    
    previousMtxSize <<- 0
    
    sampleInfo <- input$sampleDataTable_cell_clicked
    # info$value will return the Sample clicked
    
    # do nothing if not clicked yet, or the clicked cell is not in the 1st col
    #   which is (index 0 for DT object)
    if (is.null(sampleInfo$value) || sampleInfo$col != 0) return()
    
    # otherwise, change the selected tabs on the client
    updateTabsetPanel(session, "navbarPage", selected = "Variants")
    
    updateTextInput(session, "varTabInput", value = sampleInfo$value)
    
    output$varTable <- DT::renderDataTable({
      data = datatable(GetVCF(sampleInfo$value),
                       selection = list(mode = "multiple", target = "cell"),
                       rownames = FALSE) %>%
        formatStyle(columns = 2, cursor = "pointer")
    })
    
    # First renderPlot
    output$coveragePlot <- renderPlot({
      PlotCoverage(sample = sampleInfo$value)
      })
  })
  
  observeEvent(input$varTable_cells_selected, {
    # Logical to determine whether to re-paint a blank plot
    redrawBlank <- FALSE
    
    originalMatrix <- input$varTable_cells_selected

    if (!(all(is.na(originalMatrix)))) {
      # Select only values that are from the correct column in the variant 
      #   table ("positions" - col 1)
      filteredMatrix <- originalMatrix[originalMatrix[ , 2] == 1, , 
                                       drop = FALSE]
      
      # Is anything left in the filteredMatrix?
      if (!(all(is.na(filteredMatrix)))) {
        # If something is left in the filtered matrix, continue...
        # First adjust the index values in the returned matrix,
        # because the DT object is 0-indexed, but the data frame is 1-indexed
        # (add 1 to the column values in the filtered matrix)
        
        currentMtxSize <- nrow(filteredMatrix)

        # if the current size of filtered matrix is different than before, 
        #   we need to re-draw the plot with new variants
        if (previousMtxSize != currentMtxSize) {
          # Adjust the indices
          newMatrix <- cbind((filteredMatrix[ , 1]), filteredMatrix[ , 2] + 1)
          # Pull the data for selected variants from the variant call file
          vcf <- GetVCF(input$sampleDataTable_cell_clicked$value)[newMatrix]
          # Make the plot with variants identified
          output$coveragePlot <- renderPlot({
            PlotCoverage(sample = input$sampleDataTable_cell_clicked$value, 
                         positions = as.numeric(vcf))
          })
        }
        previousMtxSize <<- as.numeric(currentMtxSize)
        
      } else if (previousMtxSize != 0) {
        redrawBlank <- TRUE
        previousMtxSize <<- 0
      }
    
    } else {
      redrawBlank <- TRUE
      previousMtxSize <<- 0
    }
    
    if (redrawBlank == TRUE) {
      output$coveragePlot <- renderPlot({
        PlotCoverage(sample = input$sampleDataTable_cell_clicked$value)
      })
    }
    
  }) # end of second observeEvent
    
}

# Run the application 
shinyApp(ui = ui, server = server)



# testing how to parse out the VCF file to get info for HighlightTrack

#ht <- HighlightTrack(trackList = list(dtrack), 
                     #start = c(2000,5000), 
                     #width = c(1, 100),
                     #inBackground = FALSE)
#test <- plotTracks(list(gtrack, ht))

#testVCFFile <- GetVCF("Baldridge_1")
#startVector <- as.numeric(testVCFFile$Position) # goes to start argument
#altAlleles <- as.character(testVCFFile$Alternative) # goes to width argument
#widthVector <- map_int(.x = altAlleles, .f = nchar)
  