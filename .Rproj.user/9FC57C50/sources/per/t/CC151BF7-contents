#----- Libraries -----#

library("shiny")
library("DT")
library("Gviz")
library("data.table")
library("tidyverse")

options(shiny.sanitize.errors = FALSE)
options(ucscChromosomeNames = FALSE)

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
  
  options(useChromosomeNames = FALSE)
  
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
  
  # is positions and widths null? 
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

  
  
  # Generate a highlight track that uses info from the VCF file to mark the
  #   locations of the variants (get the width from the number of characters
  #   in the alternative allele)
  #vcfFile <- GetVCF(sample)
  #startVector <- as.numeric(vcfFile$Position) # goes to start argument
  #widthVector <- map_int(.x = as.character(vcfFile$Alternative), 
   #                      .f = nchar) # goes to width argument
  #htrack <- HighlightTrack(trackList = list(dtrack),
  #                         start = startVector,
  #                         width = widthVector,
  #                         inBackground = FALSE,
  #                         fill = "#FFE3E6b8")
  # Plot
  #coveragePlot <- plotTracks(list(gtrack, htrack))
}

#----- Set up Data Table for App -----#

# Read in sample alignment count data, calculate percent MNV
readCounts <- read.delim("../alignment_counts.txt")
readCounts <- readCounts %>%
  mutate("percent_MNV" = round(100*(primary_alignments/total_reads), 
                               digits = 2))

# Read in genome coverage count data, calculate avg. genome coverage
rawFiles <- list.files(path = "../sample_data/genome_coverage/", 
                       pattern = "*_coverage.txt")

# Some constants for calculation/formatting avg. genome coverage:
genome_size <- 7383
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
           verbatimTextOutput("cell_clicked", TRUE),
           verbatimTextOutput("cell_selected", TRUE),
           plotOutput("coveragePlot"))
  
  #tabPanel(title = "Coverage", textInput("covPlotInput", "Row ID"),
           #plotOutput("coveragePlot"))
  
)

server <- function(input, output, session) {
  
  # columnDefs - hide the 3rd column (number primary alignments) 
  #   which is index 2 in this case
  # formatStyle - add CSS style 'cursor: pointer' to the 1st column (i.e. sample)
  #   which is index 1 in this case
  output$sampleDataTable <- DT::renderDataTable({
    
    data = datatable(sampleData, selection = "single", rownames = FALSE,
                     options = list(
                       columnDefs = list(list(visible = FALSE, 
                                              targets = c(2))))) %>%
      formatStyle(columns = 1, cursor = "pointer") 
    
  })
  
  # see https://yihui.shinyapps.io/DT-click/
  
  observeEvent(input$sampleDataTable_cell_clicked, {
    
    info <- input$sampleDataTable_cell_clicked
    # info$value will return the Sample that can be passed to both GetVCF()
    #   and PlotCoverage()

    # do nothing if not clicked yet, or the clicked cell is not in the 1st col
    #   which is index 0 in this case
    if (is.null(info$value) || info$col != 0) return()
    
    # otherwise, change the selected tabs on the client
    updateTabsetPanel(session, "navbarPage", selected = "Variants")
  
    updateTextInput(session, "varTabInput", value = info$value)
    
    output$varTable <- DT::renderDataTable({
      data = datatable(GetVCF(info$value),
                       selection = list(mode = "multiple", target = "cell"),
                       rownames = FALSE) %>%
        formatStyle(columns = 2, cursor = "pointer")
    })
    
    output$coveragePlot <- renderPlot({
      withProgress(message = "drawing plot", value = 0.9, {
        PlotCoverage(sample = info$value)
      })
    })
    
    cell_selected_reactive <- eventReactive(input$varTable_cells_selected, {
      
      originalMatrix <- input$varTable_cells_selected
      
      if (!(all(is.na(originalMatrix)))) {
        
        filteredMatrix <- originalMatrix[originalMatrix[ , 2] == 1, , drop = FALSE]
        
        if (!(all(is.na(filteredMatrix)))) {
          
      
        # because of 0 - indexing
        newMatrix <- cbind((filteredMatrix[ , 1]), filteredMatrix[ , 2] + 1)
        vcf <- GetVCF(info$value)[newMatrix]
        
        output$coveragePlot <- renderPlot({
          withProgress(message = "drawing plot", value = 0.9, {
            PlotCoverage(sample = info$value, positions = as.numeric(vcf))
            })
          })
        
      } else {
        
        output$coveragePlot <- renderPlot({
          withProgress(message = "drawing plot", value = 0.9, {
            PlotCoverage(sample = info$value)
          })
        })
      }
      } else {
        
        output$coveragePlot <- renderPlot({
          withProgress(message = "drawing plot", value = 0.9, {
            PlotCoverage(sample = info$value)
          })
        })
      }
      
    })
    
    output$cell_clicked <- renderPrint(cell_selected_reactive())

  })
}

    
    
    
    # each time cell_clicked event is fired, get $value, check list for
    #   $value, add if not in list, remove if it does
    
    #output$coveragePlot <- renderPlot({
      
    #})
    
    
    
    #output$cell_clicked <- renderPrint(input$varTable_cell_clicked)
    #output$cell_selected <- renderPrint(input$varTable_cells_selected)

    #output$test <- renderPrint(input$varTable_cell_clicked)
    
    # NOW - put some kind of eventReactive or observeEvent....as long as
    #  output$test > 0, and which updates every time output$test updates
    #   and re-render the coverage plot with the highlights
    
    #observeEvent(input$varTable_cell_clicked, {
      
     # variantInfo <- input$varTable_cells_selected
      #showOutput <- eventReactive(variantInfo, {
       # info <- variantInfo
        #return(info)
      #})
      #output$test <- renderPrint(showOutput())
      # do nothing if not yet clicked, or clicked cell isn't in the 2nd col
      #   which is index 1 in this case
      #if (is.null(variantInfo$value) || variantInfo$col != 1) return()
      
      #positionSelected <- as.integer(variantInfo$value)

      # update the coverage plot with highlight tracks
      #  this needs to run different when a new sample is selected
      #altAlleleWidth <- reactive({
       # vcf <- GetVCF(info$value)
      #  vcfRow <- subset(vcf, vcf$Position == positionSelected)
      #  alternative <- vcfRow[ , "Alternative"]
      #  altWidth <- nchar(alternative)
      #  return(altWidth)
      #})
      
      #output$test <- renderPrint(variantInfo)
      
      #output$coveragePlot <- renderPlot({
      #  PlotCoverage(sample = info$value, positions = positionSelected,
      #               widths = altAlleleWidth())
      #})
      
    #}) # end of nested observeEvent()
    
  #})# end of outer-most observeEvent()
  
#}

# Run the application 
shinyApp(ui = ui, server = server)


# testing how to read in a VCF file

#sample <- "Baldridge_10"
#variantSample <- paste0("../../variants/", sample, "_variants.vcf")
#vcfFile <- read.delim(variantSample,comment.char = "#", header = FALSE, colClasses = "character")
#vcfHeaders <- c("RefGenome", "Position", "ID", "Reference", "Alternative","Quality", "Filter", "Info", "Format", "Values")
#names(vcfFile) <- vcfHeaders
# The last number in Values will be the allelic depth - unfiltered number of 
#   reads supporting the reported allele(s)
#vcfFileFormatted <- vcfFile %>%
  #mutate("AllelicDepth" = map_chr(.x = str_split(Values, pattern = ","),
   #                                .f = tail, n = 1))
#vcfFileFinal <- select(vcfFileFormatted, -c("Filter", "Info", "Format", "Values"))

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
  