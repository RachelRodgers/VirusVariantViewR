#----- Libraries -----#

library("shiny")
library("DT")
library("stringr")
library("tidyverse")

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

#----- Format, Run App -----#

ui <- navbarPage(
  
  # title id refers to navbarPage object
  title = "MNV (CR6) Variants", id = "navbarPage",
  
  tabPanel(title = "Table", DT::dataTableOutput("sampleDataTable")),
  
  tabPanel(title = "Variants", textInput("varTabInput", "Row ID"),
           DT::DTOutput("varTable"))
  
)

server <- function(input, output, session) {
  
  # columnDefs - hide the 3rd column (number primary alignments) 
  #   which is index 2 in this case
  # formatStyle - add CSS style 'cursor: pointer' to the 1st column (i.e. sample)
  #   which is index 1 in this case
  output$sampleDataTable <- DT::renderDataTable({
    
    data = datatable(sampleData, selection = "none", rownames = FALSE,
                     options = list(
                       columnDefs = list(list(visible = FALSE, 
                                              targets = c(2))))) %>%
      formatStyle(columns = 1, cursor = "pointer")
    
  })
  
  observeEvent(input$sampleDataTable_cell_clicked, {
    
    info <- input$sampleDataTable_cell_clicked
    # do nothing if not clicked yet, or the clicked cell is not in the 1st col
    #   which is index 0 in this case
    if (is.null(info$value) || info$col != 0) return()
    # otherwise, change the selected tabs on the client
    updateTabsetPanel(session, "navbarPage", selected = "Variants")
    updateTextInput(session, "varTabInput", value = info$value)
    
    output$varTable <- DT::renderDataTable({
      data = datatable(GetVCF(info$value))
    })
    
  })# end of observeEvent

}

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


  