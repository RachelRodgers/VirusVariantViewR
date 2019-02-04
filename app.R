#----- Libraries -----#

library("shiny")
library("shinythemes")
library("shinyjs")
library("DT")
library("Gviz")
library("data.table")
library("tidyverse")

source("./global.R")

options(shiny.sanitize.errors = FALSE) # need to see the error
options(ucscChromosomeNames = FALSE) # for Gvis

#--------------------- Format, Run App ---------------------#

ui <- tagList(
  
  useShinyjs(),
  
  navbarPage(
    
    theme = shinytheme("flatly"),
    title = "VariantViewR", id = "navbarpage",
    
    tabPanel(title = "Data Set Selection",
           selectInput(inputId = "dataSetSelect",
                       label = "Available Data Sets:",
                       choices = c("baldridge_rumspringa", "craig_mnv")),
           actionButton(inputId = "go", label = "Go"),
           verbatimTextOutput("buttonValue")),
    
    tabPanel(title = "Sample Table", value = "sampleTab",
           DT::dataTableOutput(outputId = "sampleDataTable"),
           actionButton(inputId = "showVariants", label = "Show Variants"),
           verbatimTextOutput("datasetValue"),
           verbatimTextOutput("sampleSelection")),
    
    tabPanel(title = "Variants", value = "variantTab",
           textInput(inputId = "varTabInput", label = "Current Sample:"),
           DT::DTOutput(outputId = "varTable"),
           plotOutput(outputId = "coveragePlot"))
  )
)


server <- function(input, output, session) {

  hideTab(inputId = "navbarpage", target = "sampleTab")
  hideTab(inputId = "navbarpage", target = "variantTab")
  
  #----- Data Set Selection -----#
  
  # Once action button is pressed, show the sample table tab
  observeEvent(input$go, {
    # move to sample table tab:
    showTab(inputId = "navbarpage", target = "sampleTab", select = TRUE)
    hideTab(inputId = "navbarpage", target = "variantTab")
    # show the name of the currently selected data set:
    #output$datasetValue <- renderText(input$dataSetSelect)
    
    # get the correct data set information to display:
    sampleData <- GenerateAlignmentCounts(dataSet = input$dataSetSelect)
    
    # columnDefs - hide the 3rd column (number primary alignments) 
    #   which is index 2 in the DT object
    # formatStyle - add CSS style 'cursor: pointer' to the 1st column (i.e. sample)
    #   which is index 1 in in the DT object
    output$sampleDataTable <- DT::renderDataTable({
      data = datatable(sampleData, 
                       selection = list(mode = "single", target = "cell"), 
                       rownames = FALSE,
                       options = list(
                         columnDefs = list(list(visible = FALSE, 
                                                targets = c(2))))) %>%
        formatStyle(columns = 1, cursor = "pointer") 
    })
  }) # end of data set selection
  
  #----- Sample Selection -----#
  
  observeEvent(input$sampleDataTable_cell_clicked, {
    # what's selected?
    sampleInfo <- input$sampleDataTable_cell_clicked
    output$datasetValue <- renderPrint(sampleInfo)
    # Do nothing if nothing has been clicked, or the clicked cell isn't in the
    #   first column (which is index 0 for DT objects)
    observeEvent(input$sampleDataTable_cells_selected, {
      cellsSelected <- input$sampleDataTable_cells_selected
      output$sampleSelection <- renderPrint(cellsSelected)
      if (all(is.na(input$sampleDataTable_cells_selected)) || 
          is.null(sampleInfo$value) || sampleInfo$col != 0) {
        # disable the button and re-hide the variant tab
        shinyjs::disable("showVariants")
        hideTab(inputId = "navbarpage", target = "variantTab")
        return()
      } else {
        shinyjs::enable("showVariants")
      }
    })
    
  })
    
  # click the button:
  observeEvent(input$showVariants, {
    # make the variant tab visible & switch to it
    showTab(inputId = "navbarpage", target = "variantTab", select = TRUE)
    sampleInfo <- input$sampleDataTable_cell_clicked
    updateTextInput(session, inputId = "varTabInput", value = sampleInfo$value)
    
    # render the variant table:
    output$varTable <- DT::renderDataTable({
      data = datatable(GetVCF(dataSet = input$dataSetSelect,
                              sample = sampleInfo$value),
                       selection = list(mode = "multiple", target = "cell"),
                       options = list(pageLength = 5),
                       rownames = FALSE) %>%
        formatStyle(columns = 2, cursor = "pointer")
    })
    
    # First renderPlot (no tracks)
    output$coveragePlot <- renderPlot({
      PlotCoverage(dataSet = input$dataSetSelect, sample = sampleInfo$value)
    })

    # global for determining when to re-draw coverage plot
    previousMtxSize <<- 0
  }) # end of sample selection observeEvent

  #----- Variant Selection -----#
  
  observeEvent(input$varTable_cells_selected, {
    # Logical to determine whether to re-paint a blank plot
    redrawBlank <- FALSE
    # What's selected on the variant table?
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
          vcf <- GetVCF(dataSet = input$dataSetSelect,
                        sample = input$sampleDataTable_cell_clicked$value)[newMatrix]
          # Make the plot with variants identified
          output$coveragePlot <- renderPlot({
            PlotCoverage(dataSet = input$dataSetSelect,
                         sample = input$sampleDataTable_cell_clicked$value, 
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
        PlotCoverage(dataSet = input$dataSetSelect,
                     sample = input$sampleDataTable_cell_clicked$value)
      })
    }
  }) # end of variant selection observeEvent
  
}

shinyApp(ui = ui, server = server)
