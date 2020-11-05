packages<-c('shiny','shinythemes','shinydashboard','shinycssloaders','readxl','xlsx','ggplot2','rapportools','plotly','dplyr')
for (package in packages){
  if(package %in% rownames(installed.packages()) == FALSE) {
    install.packages(package)}
}
library(shiny)
library(shinythemes)
library(shinydashboard)
library(readxl)
library(xlsx)
library(ggplot2)
library(rapportools)
library(plotly)
library(dplyr)


source('global.R', local = TRUE)
ui<- dashboardPage(
  dashboardHeader(title = 'Bowl\'s project'),
  dashboardSidebar(
    sidebarMenu(
      menuItem('Data', tabName = 'start', icon = icon('chart-line')),
      menuItem("How to use", tabName = "guideline", icon = icon('map'))
    ),
    fileInput(inputId = 'file',label = 'Select input file:',multiple = FALSE),
    tableOutput('filename'),
    tags$hr(),
    actionButton(inputId = "analyze", label = "Get data")
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "guideline",
              navbarPage(title = 'About this website'),
              tags$div(
                tags$p('This website is customerized for WanYing\'s realtime PCR results analysis. Please make sure you have a correct format before you upload your file.')
              ),
              tags$h3('1. Check your file format'),
              tags$div(
                tags$p('Store your data in the sheet named, "Results".'),
                tags$br(),
                tags$img(src='excel_format.png'),
                tags$br(),
                tags$p('Make sure your data starts from row 38'),
                tags$p('Row 38th will be the column names and only data in column "Well Position"
                        ,"Sample Name","Target Name",8,9,10,15,23 will be captured')
              ),
              tags$h3('2. Select and upload your file on the Sidebar')
      ),
      tabItem(tabName = 'start',
              navbarPage(title ='Data analysis',
                         tabPanel('Set up Ctrls', icon = icon('calendar-plus'),
                                  box(width = 4,
                                      uiOutput("controls")
                                  ),
                                  box(width = 4,
                                      checkboxGroupInput("selected_controls", "Selected Control"),
                                      uiOutput('internal_control'),
                                      uiOutput("show_results"),
                                      actionButton(inputId = "analyze_ctrl", label = "Analyze")
                                  ),
                                  box(width = 4, title = 'Control data:',
                                      tableOutput('ctrls'),
                                      downloadLink("download", "Download xlsx file")
                                  )
                         ),
                         tabPanel('Raw data', icon = icon('file'),
                                  tableOutput('rawdata')
                         ),
                         tabPanel('filtered', icon = icon('calendar'),
                                  tags$p('In this table, some ineligible data was filtered out whose Ct.SD > 0.1 and the ratio of CT among individual samples is greater than 1.1 or lesser than 0.9.'),
                                  tags$br(),
                                  tableOutput('filtered')
                         ),
                         tabPanel('Calibration', icon = icon('calendar'),
                                  tableOutput('calibration')
                         ),
                         tabPanel('Relative Data', icon = icon('calendar-check'),
                                  uiOutput('tables'),
                         ),
                         tabPanel('Plots', icon = icon('chart-bar'),
                                  radioButtons(inputId = 'plottype',label = 'Please select output plot type:',choices = c('Dot-Bar-plot','Dot-plot', 'Bar-plot','Box-plot','Dot-Box-plot' ), selected = 'Dot-Bar-plot'),
                                  uiOutput("plotly"),
                         ),
                         tabPanel('Summary', icon = icon('signal'),
                                  uiOutput('summary')
                                  
                         )
              )
      )
      
    )
  ),
  
)


server <- function(input, output, session){
  source('global.R', local = TRUE)
  # functions and variables
  readfile<- eventReactive(input$analyze,{
    read_data(input$file$datapath)
  })
  filtered_data<-eventReactive(input$analyze,{
    get_filtered_data(readfile())
  })
  all_sample_data<-eventReactive(input$analyze,{
    get_all_sample_data(filtered_data())
  })
  relative_data_intra<-eventReactive(input$analyze_ctrl,{
    get_relative_samples_data(all_sample_data())
  })
  
  Ctrls.Mean.dCt_table<-eventReactive(input$analyze_ctrl,{
    get_Ctrls.Mean.dCt_table(relative_data_intra())
  })
  
  relative_data_inter<-eventReactive(input$analyze_ctrl,{
    get_relative_data_inter(relative_data_intra())
  })

  # Output
  # Create choosable control list
  output$controls <- renderUI({
    groups <-row.names(table(readfile()$Sample.Name))
    Ctrl.Samples<-set_control_groups(groups)
    checkboxGroupInput("controls", "Choose Controls:", groups, selected = Ctrl.Samples)
  })
  output$show_results <- renderUI({
    groups <-row.names(table(relative_data_inter()$Target.Name))
    checkboxGroupInput("show_results", "Look up genes in groups:", groups, selected = groups)
  })
  output$internal_control<-renderUI({
    control <-detect_internal_gene(readfile())
    selectInput('internal_control', 'Please select internal control', choices = c('GAPDH','ACTIN','HPRT'), selected = control)
  })
  # Session
  observe({
    x <- input$controls
    
    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)
    
    # Can also set the label and select items
    updateCheckboxGroupInput(session, "selected_controls",
                             label = paste("Controls are", length(x)),
                             choices = x,
                             selected = x
    )
  })
  output$filename<-renderTable({
    input$file[1]
  })
  output$rawdata<-renderTable({
    readfile()
  })
  output$filtered<-renderTable({
    filtered_data()
  })  
  output$calibration<-renderTable({
    all_sample_data()
  })
  output$ctrls<-renderTable({
    Ctrls.Mean.dCt_table()
  })

  output$tables<-renderUI({
    table_output_list <- lapply(1:length(input$show_results), function(i) {
      tablename <- paste0("table", i)
      table_output_object <- plotlyOutput(tablename)
      table_output_object <- renderTable({
        source('global.R', local = TRUE)
        gene_expression_table(input$show_results[i]) 
      })
    })
    return(table_output_list)
  })
  output$plotly<-renderUI({
    plot_output_list <- lapply(1:length(input$show_results), function(i) {
      plotname <- paste0("plotly", i)
      plot_output_object <- plotlyOutput(plotname)
      plot_output_object <- renderPlotly({
        source('global.R', local = TRUE)
        draw_plot(input$show_results[i],input$plottype) 
      })
    })
    return(plot_output_list)
  })
  output$summary<-renderUI({
    table_output_list <- lapply(1:length(input$show_results), function(i) {
      tablename <- paste0("table", i)
      table_output_object <- plotlyOutput(tablename)
      table_output_object <- renderTable({
      target.gene<-group_by(relative_data_inter()[relative_data_inter()$Target.Name==input$show_results[i],], Group)
      df<-summarize(target.gene, count = n(), 
              Max=max(exprs) , Min=min(exprs),
              Q1=quantile(exprs, na.rm = T)[2], Q2=quantile(exprs, na.rm = T)[3], Q3=quantile(exprs, na.rm = T)[4],
              Mean = mean(exprs, na.rm = T), SD = sd(exprs, na.rm = T))
      names(df)<-c(input$show_results[i], 'Number','Max','Min','Q1','Q2','Q3', 'Mean','SD')
      return(df)
      })
    })
    return(table_output_list)
  })
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$file[1], ".summary.xlsx")
    },
    content = function(file) {
        write.xlsx(readfile(), sheetName = 'raw data',file)
        write.xlsx(filtered_data(), sheetName = 'filtered',file)
        write.xlsx(Ctrls.Mean.dCt_table(), sheetName = 'Control',file, append = TRUE)
        write.xlsx(all_sample_data(), sheetName = 'Calibration',file, append = TRUE)
        write.xlsx(relative_data_inter(), sheetName = 'Relative exprs',file, append = TRUE)
    })
}

shinyApp(ui, server)
