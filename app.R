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
                       choices = c("craig_mnv", "baldridge_rumspringa")),
           actionButton(inputId = "go", label = "Go"),
           verbatimTextOutput("buttonValue")),
    
    tabPanel(title = "Sample Table", value = "sampleTab",
           DT::dataTableOutput(outputId = "sampleDataTable"),
           actionButton(inputId = "showVariants", label = "Show Variants"),
           verbatimTextOutput("datasetValue"),
           verbatimTextOutput("sampleSelection"),
           verbatimTextOutput("filteredMatrix")),
    
    tabPanel(title = "Variants", value = "variantTab",
           textInput(inputId = "varTabInput", label = "Current Sample:"),
           uiOutput(outputId = "variantTables"),
           uiOutput(outputId = "coveragePlots"))#,
           #div(id = "placeholder"))#,
           #DT::DTOutput(outputId = "varTable"),
           #plotOutput(outputId = "coveragePlot"))
  )
)


server <- function(input, output, session) {

  hideTab(inputId = "navbarpage", target = "sampleTab")
  sampleVec <- vector(mode = "character")
  hideTab(inputId = "navbarpage", target = "variantTab")
  
  #----- Data Set Selection -----#
  
  # Once action button is pressed, show the sample table tab
  observeEvent(input$go, {
    # move to sample table tab:
    showTab(inputId = "navbarpage", target = "sampleTab", select = TRUE)
    sampleVec <<- vector(mode = "character")
    hideTab(inputId = "navbarpage", target = "variantTab")
    # show the name of the currently selected data set:
    #output$datasetValue <- renderText(input$dataSetSelect)
    
    # get the correct data set information to display:
    sampleData <- GenerateAlignmentCounts(dataSet = input$dataSetSelect)
    
    # columnDefs - hide the 3rd column (number primary alignments) 
    #   which is index 2 in the DT object
    # formatStyle - add CSS style 'cursor: pointer' to the 1st column (i.e. sample)
    #   which is index 1 in in the DT object
    # Table set up to allow multiple cell selections
    output$sampleDataTable <- DT::renderDataTable({
      data = datatable(sampleData, 
                       selection = list(mode = "multiple", target = "cell"), 
                       rownames = FALSE,
                       options = list(
                         columnDefs = list(list(visible = FALSE, 
                                                targets = c(2))))) %>%
        formatStyle(columns = 1, cursor = "pointer") 
    })
  }) # end of data set selection
  
  #----- Sample Selection -----#
  
  # Select one or more samples, store the sample names in a chr vector.

  GetSampleVec <- eventReactive(input$sampleDataTable_cell_clicked, {
    currentCell <- input$sampleDataTable_cell_clicked
    currentCellCol <- currentCell$col
    if (currentCellCol == 0) {
      currentSample <- currentCell$value
      if (!(currentSample %in% sampleVec)) {
        sampleVec <<- c(sampleVec, currentSample)
      } else { # remove sample
        sampleVec <<- sampleVec[sampleVec != currentSample]
      }
    }
    return(sampleVec)
  }, ignoreNULL = FALSE)
  
  # show value of sample vector in the sample table tab
  #observe({
   # output$filteredMatrix <- renderPrint(GetSampleVec())
    #})
  
  # determine when to enable/disable the "Show Variants" button
  SelectedSampleIndices <- eventReactive(input$sampleDataTable_cell_clicked, {
    sampleMtxOriginal <- input$sampleDataTable_cells_selected
    sampleMtxFiltered <- NULL
    # if sampleMtxOriginal is not empty, filter out incorrect cells
    if (!(all(is.na(sampleMtxOriginal)))) {
      sampleMtxFiltered <- sampleMtxOriginal[sampleMtxOriginal[, 2] == 0, ,
                                             drop = FALSE]
    }
    return(sampleMtxFiltered)
  })
  
  # enables/disables the Show Variants button:
  observeEvent(input$sampleDataTable_cell_clicked, {
    # Do nothing if nothing has been clicked, or the clicked cell isn't in the
    #   first column (which is index 0 for DT objects)
    if (is.null(SelectedSampleIndices()) | all(is.na(SelectedSampleIndices()))) {
      shinyjs::disable("showVariants")
      hideTab(inputId = "navbarPage", target = "variantTab")
      return()
      } else {
          shinyjs::enable("showVariants")
        }
  })
  
  observeEvent(input$sampleDataTable_cell_clicked, {
    # what's selected?
    # last clicked cell:
    output$datasetValue <- renderPrint(input$sampleDataTable_cell_clicked)
    # current filtered Matrix
    output$sampleSelection <- renderPrint(SelectedSampleIndices())
    # current filtered Matrix
    #output$filteredMatrix <- renderPrint(GetSampleVec())
    
  })
  
  
  #----- View Variants -----#  
  
  # click the button:
  observeEvent(input$showVariants, {
    # make the variant tab visible & switch to it
    showTab(inputId = "navbarpage", target = "variantTab", select = TRUE)
    #sampleInfo <- input$sampleDataTable_cell_clicked
    #updateTextInput(session, inputId = "varTabInput", value = sampleInfo$value)
    updateTextInput(session, inputId = "varTabInput", value = GetSampleVec())
    
    output$variantTables <- renderUI({
      allVariantData <- map(isolate(GetSampleVec()), GetVCF, 
                            dataSet = input$dataSetSelect)
      allVariantTables <- map(allVariantData, function(x) {
        DT::renderDataTable(expr = 
                              (data = datatable(x,
                                                selection = list(mode = "multiple", 
                                                                 target = "cell"),
                                                options = list(pageLength = 5), 
                                                rownames = FALSE)) %>%
                              formatStyle(columns = 2, cursor = "pointer"))
        })
    })
    
    output$coveragePlots <- renderUI({
      allCoveragData <- map(GetSampleVec(), function(y) {
        renderPlot({PlotCoverage(dataSet = input$dataSetSelect, sample = y)})
      })

  })
  })
  
    
    # First renderPlot (no tracks)
    #output$coveragePlot <- renderPlot({
      #PlotCoverage(dataSet = input$dataSetSelect, sample = sampleInfo$value)
    #})

    # global for determining when to re-draw coverage plot
    #previousMtxSize <<- 0
  #}) # end of sample selection observeEvent

  
}

shinyApp(ui = ui, server = server)



# test loop
#testVector <- c("Baldridge_1", "Baldridge_2")
#listOfTracks <- map(testVector, PlotCoverage, dataSet = "craig_mnv")

#testFunction <- function(track_list) {
#  for (i in 1:length(track_list)) {
#    plotTracks(trackList = list(track_list[[i]][[1]],
#                                track_list[[i]][[2]]))
#  }
#}

#plotTest <- testFunction(listOfTracks)

