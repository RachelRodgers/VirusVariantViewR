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
           actionButton(inputId = "showCoverage", label = "Show Coverage"),
           verbatimTextOutput("verbatimOutput1"),
           verbatimTextOutput("verbatimOutput2"),
           verbatimTextOutput("verbatimOutput3")),
    
    tabPanel(title = "Coverage", value = "coverageTab",
           textInput(inputId = "varTabInput", label = "Current Sample:"),
           h4("Display Options:"),
           checkboxInput(inputId = "showVarTables", 
                         label = "Show Variant Tables",
                         value = FALSE),
           checkboxInput(inputId = "showCovPlots",
                         label = "Show Coverage Plots",
                         value = TRUE),
           uiOutput(outputId = "variantTables"),
           uiOutput(outputId = "coveragePlots")),
    
    tabPanel(title = "Variant Inspector", value = "variantInsp",
             textInput(inputId = "variantInspTextBox", label = "Current Sample:"),
             DT::DTOutput(outputId = "variantInspDT"),
             plotOutput(outputId = "variantInspCovPlot"))
  )
)


server <- function(input, output, session) {

  hideTab(inputId = "navbarpage", target = "sampleTab")
  hideTab(inputId = "navbarpage", target = "coverageTab")
  #hideTab(inputId = "navbarpage", target = "variantInsp")
  
  # define the vector that will hold user-selected samples for coverage display:
  sampleVec <- vector(mode = "character") 
  
  #----- Data Set Selection -----#
  
  observeEvent(input$go, {
    showTab(inputId = "navbarpage", target = "sampleTab", select = TRUE)
    hideTab(inputId = "navbarpage", target = "coverageTab")
    
    # clear the vector that will hold user-selected samples for coverage disaplay:
    sampleVec <<- vector(mode = "character")
    
    #----- Display Sample Data Table -----#
    
    sampleData <- GenerateAlignmentCounts(dataSet = input$dataSetSelect)
    
    # columnDefs - hide the 3rd column (index 2) (number primary alignments)
    # formatStyle - add CSS style 'cursor: pointer' to the 1st column (i.e. sample)
    #   which is index 1 in in the DT object
    # Table set up to allow multiple cell selections:
    output$sampleDataTable <- DT::renderDataTable({
      data = datatable(sampleData, 
                       selection = list(mode = "multiple", target = "cell"), 
                       rownames = FALSE,
                       options = list(
                         columnDefs = list(list(visible = FALSE, 
                                                targets = c(2))),
                         pageLength = 10,
                         bLengthChange = 0)) %>%
        formatStyle(columns = 1, cursor = "pointer") 
    })
  })
  
  #----- Sample Selection -----#

  GetSampleVec <- eventReactive(input$sampleDataTable_cells_selected, {
    # Returns character vector of Sample IDs.
    currentCell <- input$sampleDataTable_cell_clicked # row, col, value info
    currentCellCol <- currentCell$col
    if (currentCellCol == 0) { # value from Sample column must be selected
      currentSample <- currentCell$value
      if (!(currentSample %in% sampleVec)) {
        sampleVec <<- c(sampleVec, currentSample)
      } else { # remove sample
        sampleVec <<- sampleVec[sampleVec != currentSample]
      }
    }
    return(sampleVec)
  })
  
  # show value of sample vector in the sample table tab
  observeEvent(input$sampleDataTable_cells_selected, {
    output$verbatimOutput1 <- renderPrint(GetSampleVec())
    })
  
  # Calculate a matrix of correctly-selected cell indices.
  #   Used to determine when to enable/disable the "Show Variants" button by seeing if
  #   anything is selected on the sample data table.  I have to do this because
  #   I can't figure out how to prevent crashing from GetSampleVec()'s initial 
  #   empty/NULL value.  Solving that would simplify the code a little bit.
  SelectedSampleIndices <- eventReactive(input$sampleDataTable_cells_selected, {
    # User is free to select "incorrect" cells so they need to be filtered first.
    sampleMtxOriginal <- input$sampleDataTable_cells_selected
    sampleMtxFiltered <- NULL
    if (!(all(is.na(sampleMtxOriginal)))) {
      sampleMtxFiltered <- sampleMtxOriginal[sampleMtxOriginal[, 2] == 0, ,
                                             drop = FALSE]
    }
    return(sampleMtxFiltered)
  })
  
  # enables/disable the Show Variants button:
  observeEvent(input$sampleDataTable_cells_selected, {
    # Do nothing if nothing has been clicked, or the clicked cell isn't in the
    #   first column (which is index 0 for DT objects)
    if (is.null(SelectedSampleIndices()) | all(is.na(SelectedSampleIndices()))) {
      shinyjs::disable("showVariants")
      hideTab(inputId = "navbarPage", target = "coverageTab")
      return()
    } else {
      shinyjs::enable("showVariants")
    }
  })
  
  observe({
    # what's selected?
    # indices of ALL selected cells
    output$verbatimOutput2 <- renderPrint(input$sampleDataTable_cells_selected)
    # indices of CORRECT selected cells
    output$verbatimOutput3 <- renderPrint(SelectedSampleIndices())
  })
  
  
  #----- View Coverage Plots & Variant Tables -----#  
  
  observeEvent(input$showCoverage, {
    # make the variant tab visible & switch to it
    showTab(inputId = "navbarpage", target = "coverageTab", select = TRUE)
    # reset the checkboxes - don't show variant tables, do show cov plots
    updateCheckboxInput(session, inputId = "showVarTables", value = FALSE)
    updateCheckboxInput(session, inputId = "showCovPlots", value = TRUE)

    updateTextInput(session, inputId = "varTabInput", value = GetSampleVec())
    
    # for showing/hiding variant tables
    observeEvent(input$showVarTables, {
      if(input$showVarTables == TRUE) {
        shinyjs::show(id = "variantTables")
      } else {
        shinyjs::hide(id = "variantTables")
      }
    })
    
    # for showing/hiding coverage plots
    observeEvent(input$showCovPlots, {
      if(input$showCovPlots == TRUE) {
        shinyjs::show(id = "coveragePlots")
      } else {
        shinyjs::hide(id = "coveragePlots")
      }
    })
    
    # render the variant tables and coverage plots dynamically:
    output$variantTables <- renderUI({
      # Variant Tables:
      allVariantData <- map(isolate(GetSampleVec()), GetVCF, 
                            dataSet = input$dataSetSelect)
      
      allVariantTables <- map(allVariantData, function(x) {
        DT::renderDataTable(expr =
                              (data = datatable(x,
                                                caption = htmltools::tags$caption(
                                                  style = 'caption-side: top; text-align: center; color:black;',
                                                  htmltools::strong(comment(x))),
                                                selection = list(mode = "multiple",
                                                                 target = "cell"),
                                                options = list(paging = FALSE,
                                                               #pageLength = 5,
                                                               bLengthChange = 0,
                                                               bFilter = 0),
                                                rownames = FALSE)) %>%
                              formatStyle(columns = 2, cursor = "pointer"))
        })
    })
    
    # Variant Plots:
    output$coveragePlots <- renderUI({
      allCoveragData <- map(GetSampleVec(), function(y) {
        renderPlot({PlotCoverage(dataSet = input$dataSetSelect, sample = y)},
                   height = 300)
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
