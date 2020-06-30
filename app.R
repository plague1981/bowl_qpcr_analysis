library(shiny)
library(shinythemes)
library(shinydashboard)
library(readxl)
library(xlsx)
library(ggplot2)
library(rapportools)
library(plotly)

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
                         tabPanel('Raw data', icon = icon('file'),
                                  tableOutput('rawdata')
                         ),
                         tabPanel('Calibration', icon = icon('calendar'),
                                  tableOutput('calibration')
                         ),
                         tabPanel('Set up Ctrls', icon = icon('calendar-plus'),
                                  box(width = 4,
                                      uiOutput("controls")
                                  ),
                                  box(width = 4,
                                      checkboxGroupInput("selected_controls", "Selected Control"),
                                      selectInput('internal_control', 'Please select internal control', choices = c('GAPDH','ACTIN','HPRT'), selected = 'GAPDH'),
                                      actionButton(inputId = "analyze_ctrl", label = "Analyze")
                                  ),
                                  box(width = 4,
                                      tableOutput('ctrls')
                                  )
                         ),
                         tabPanel('Relative Data', icon = icon('calendar-check'),
                                  uiOutput("show_results"),
                                  tableOutput('relative_inter_test'),
                                  tableOutput('relative_inter_1'),
                                  tags$hr(),
                                  tableOutput('relative_inter_2'),
                                  tags$hr(),
                                  tableOutput('relative_inter_3'),
                                  tags$hr(),
                                  tableOutput('relative_inter_4'),
                                  tags$hr(),
                                  tableOutput('relative_inter_5'),
                                  tags$hr(),
                                  tableOutput('relative_inter_6'),
                         ),
                         tabPanel('Plots', icon = icon('chart-bar'),
                                  radioButtons(inputId = 'plottype',label = 'Please select output plot type:',choices = c('Dot-plot', 'Box-plot','Dot-Box-plot'),selected = 'Dot-plot'),
                                  plotOutput('relative_inter_1_plot'),
                                  tags$hr(),
                                  plotOutput('relative_inter_2_plot'),
                                  tags$hr(),
                                  plotOutput('relative_inter_3_plot'),
                                  tags$hr(),
                                  plotOutput('relative_inter_4_plot'),
                                  tags$hr(),
                                  plotOutput('relative_inter_5_plot'),
                                  tags$hr(),
                                  plotOutput('relative_inter_6_plot'),
                                  tags$hr(),
                         ),
                         tabPanel('Statistic', icon = icon('signal'),
                                  uiOutput("plots")
                                  )
              )
      )
      
    )
  ),
  
)


server <- function(input, output, session){
  # functions and variables
  readfile<- eventReactive(input$analyze,{ 
    file2table<-read_excel(input$file$datapath, col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[-(1:37),]
    colnames(file2table)<-read_excel(input$file$datapath, col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[37,]
    file2table<-data.frame(file2table)[,c(1,3,4,8,9,10,15,23)]
    return(file2table)
  })
  all_sample_data<-eventReactive(input$analyze,{ 
    all_sample_data<-NULL
    for (sample_name in row.names(table(readfile()[,'Sample.Name']))){
      sample_A<-readfile()[readfile()$Sample.Name==sample_name,]
      gene_data<-NULL
      for (gene in rownames(table(sample_A$Target.Name))){
        gene_data<-rbind(gene_data, sample_A[sample_A$Target.Name==gene, ][1,c(2,3,5)])
      }
      
      new_Ct.Mean_list<-NULL
      lq_gene<-NULL
      for (gene in rownames(table(sample_A$Target.Name))){
        if (sample_A[sample_A$Target.Name==gene, ]$Ct.Threshold[1]>0.1){
          x<-as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[1])
          y<-as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[2])
          z<-as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[3])
          a<-abs(x-y)
          b<-abs(y-z)
          c<-abs(x-z)
          if (min(c(a,b,c))==a){
            new_Ct.Mean<-(as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[1])+as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[2]))/2
            new_Ct.Mean_list<-c(new_Ct.Mean_list,new_Ct.Mean)
          } else if (min(c(a,b,c))==b){
            new_Ct.Mean<-(as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[2])+as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[3]))/2
            new_Ct.Mean_list<-c(new_Ct.Mean_list,new_Ct.Mean)
          } else if (min(c(a,b,c))==c){
            new_Ct.Mean<-(as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[1])+as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[3]))/2
            new_Ct.Mean_list<-c(new_Ct.Mean_list,new_Ct.Mean)
          }
          lq_gene<-c(lq_gene,gene)
        } 
      }
      df<-NULL
      df$new_Ct_table<-data.frame(lq_gene)
      df<-cbind(df, data.frame(new_Ct.Mean_list))
      
      for (n in 1:nrow(df)){
        for (m in 1:nrow(gene_data)){
          if (gene_data[m,'Target.Name']==df[n,'lq_gene']){
            gene_data[m,'Ct.Mean']<-df[n,'new_Ct.Mean_list']
          }
        }
      }
      all_sample_data<-rbind(all_sample_data,gene_data)
    }
    return(all_sample_data)
  })
  relative_data_intra<-eventReactive(input$analyze_ctrl,{
    # use relative_samples_data.R
    #source('relative_samples_data.R', local = TRUE)
    relative_samples_data<-NULL
    for (sample_name in row.names(table(all_sample_data()$Sample.Name))){
      A_sample_data<-all_sample_data()[all_sample_data()$Sample.Name==sample_name,]
      Crtl_Ct<-A_sample_data[A_sample_data$Target.Name==input$internal_control,]$Ct.Mean
      genes<-NULL
      genes_dCt<-NULL
      relative_sample_data<-NULL
      for (gene in A_sample_data$Target.Name){
        if (gene!=input$internal_control){
          Gene_Ct<-A_sample_data[A_sample_data$Target.Name==gene,]$Ct.Mean
          gene_dCt<-(-(as.numeric(Gene_Ct)-as.numeric(Crtl_Ct)))
          genes<-c(genes,gene)
          genes_dCt<-c(genes_dCt,gene_dCt)
        }
        relative_sample_data<-data.frame(rep(sample_name,length(genes)))
        relative_sample_data<-cbind(relative_sample_data, genes, genes_dCt)
      }
      relative_samples_data<-rbind(relative_samples_data,relative_sample_data)
    }
    # add column names
    colnames(relative_samples_data)<-c('Sample.Name','Target.Name','dCt')
    return(relative_samples_data)
  })
  
  Ctrls.Mean.dCt_table<-eventReactive(input$analyze_ctrl,{
    
    Samples.Name<-row.names(table(relative_data_intra()$Sample.Name))
    genes<-row.names(table(relative_data_intra()$Target.Name))
    Ctrl.Samples<-input$controls
    Ctrls.Mean.dCt<-NULL
    for (gene in genes ) {
      Ctrls.dCt<-NULL
      for (Ctrl.Sample in Ctrl.Samples){
        Ctrl.dCt<-relative_data_intra()[relative_data_intra()$Sample.Name==Ctrl.Sample,][relative_data_intra()[relative_data_intra()$Sample.Name==Ctrl.Sample,]$Target.Name==gene,]$dCt
        Ctrls.dCt<-c(Ctrls.dCt,Ctrl.dCt)
      }
      Ctrl.Mean.dCt<-mean(Ctrls.dCt)
      Ctrls.Mean.dCt<-c(Ctrls.Mean.dCt,Ctrl.Mean.dCt)
    }
    Ctrls.Mean.dCt_table<-data.frame(genes)
    Ctrls.Mean.dCt_table<-cbind(Ctrls.Mean.dCt_table,Ctrls.Mean.dCt)
    colnames(Ctrls.Mean.dCt_table)<-c('Target.Name','dCt')
    return(Ctrls.Mean.dCt_table)
  })
  
  relative_data_inter<-eventReactive(input$analyze_ctrl,{
    relative_data_inter<-relative_data_intra()
    ddCt_list<-NULL
    for (n in 1:nrow(relative_data_inter)){
      if (relative_data_inter[n,]$Target.Name=='BSP'){
        ddCt<-(relative_data_inter[n,'dCt']-Ctrls.Mean.dCt_table()[Ctrls.Mean.dCt_table()$Target.Name=='BSP',]$dCt)
      } else if (relative_data_inter[n,]$Target.Name=='Col1a1') {
        ddCt<-(relative_data_inter[n,'dCt']-Ctrls.Mean.dCt_table()[Ctrls.Mean.dCt_table()$Target.Name=='Col1a1',]$dCt)
      } else if (relative_data_inter[n,]$Target.Name=='Col2a1') {
        ddCt<-(relative_data_inter[n,'dCt']-Ctrls.Mean.dCt_table()[Ctrls.Mean.dCt_table()$Target.Name=='Col2a1',]$dCt)
      } else if (relative_data_inter[n,]$Target.Name=='OSX') {
        ddCt<-(relative_data_inter[n,'dCt']-Ctrls.Mean.dCt_table()[Ctrls.Mean.dCt_table()$Target.Name=='OSX',]$dCt)
      } else if (relative_data_inter[n,]$Target.Name=='Runx2') {
        ddCt<-(relative_data_inter[n,'dCt']-Ctrls.Mean.dCt_table()[Ctrls.Mean.dCt_table()$Target.Name=='Runx2',]$dCt)
      } else if (relative_data_inter[n,]$Target.Name=='SOX9') {
        ddCt<-(relative_data_inter[n,'dCt']-Ctrls.Mean.dCt_table()[Ctrls.Mean.dCt_table()$Target.Name=='SOX9',]$dCt)
      } else if (relative_data_inter[n,]$Target.Name=='Col10a1'){
        ddCt<-(relative_data_inter[n,'dCt']-Ctrls.Mean.dCt_table()[Ctrls.Mean.dCt_table()$Target.Name=='Col10a1',]$dCt)
      }
      ddCt_list<-c(ddCt_list,ddCt)
    }
    
    exprs<-2^ddCt_list
    relative_data_inter$ddCt<-ddCt_list
    relative_data_inter$exprs<-exprs
    group_list<-NULL
    
    for (n in 1:nrow(relative_data_inter)){
      if (grepl('^NC',relative_data_inter$Sample.Name[n])& grepl('D0$',relative_data_inter$Sample.Name[n])){
        group<-'NC D0'
      } else if (grepl('^NC',relative_data_inter$Sample.Name[n])& grepl('D14$',relative_data_inter$Sample.Name[n])){
        group<-'NC D14'
      } else if (grepl('^NC',relative_data_inter$Sample.Name[n])& grepl('D21$',relative_data_inter$Sample.Name[n])){
        group<-'NC D21'
      } else if (grepl('^AB',relative_data_inter$Sample.Name[n])& grepl('D0$',relative_data_inter$Sample.Name[n])){
        group<-'AB D0'
      } else if (grepl('^AB',relative_data_inter$Sample.Name[n])& grepl('D14$',relative_data_inter$Sample.Name[n])){
        group<-'AB D14'
      } else if (grepl('^AB',relative_data_inter$Sample.Name[n])& grepl('D21$',relative_data_inter$Sample.Name[n])){
        group<-'AB D21'
      } else if (grepl('^wt',relative_data_inter$Sample.Name[n])& grepl('Undiffer$',relative_data_inter$Sample.Name[n])){
        group<-'WT Undiff'
      } else if (grepl('^wt',relative_data_inter$Sample.Name[n])& grepl('Differ$',relative_data_inter$Sample.Name[n])){
        group<-'WT Diff'
      } else if (grepl('^mu',relative_data_inter$Sample.Name[n])& grepl('Undiffer$',relative_data_inter$Sample.Name[n])){
        group<-'MU Undiff'
      } else if (grepl('^mu',relative_data_inter$Sample.Name[n])& grepl('Differ$',relative_data_inter$Sample.Name[n])){
        group<-'MU Diff'
      }
      group_list<-c(group_list, group)
    }
    
    relative_data_inter$Group<-group_list
    return(relative_data_inter[order(relative_data_inter$Target.Name,relative_data_inter$Group,decreasing = TRUE),])
  })
  # Output
  # Create choosable control list
  output$controls <- renderUI({
    groups <-row.names(table(readfile()$Sample.Name))
    checkboxGroupInput("controls", "Choose Controls:", groups)
  })
  output$show_results <- renderUI({
    groups <-row.names(table(relative_data_inter()$Target.Name))
    checkboxGroupInput("show_results", "Look up genes in groups:", groups, selected = groups)
  })
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
  output$calibration<-renderTable({
    all_sample_data()
  })
  output$ctrls<-renderTable({
    Ctrls.Mean.dCt_table()
  })
  output$relative_inter_1<-renderTable({
    source('global.R', local = TRUE)
    gene_expression_table(input$show_results[1])
  })
  output$relative_inter_2<-renderTable({
    source('global.R', local = TRUE)
    gene_expression_table(input$show_results[2])
  })
  output$relative_inter_3<-renderTable({
    source('global.R', local = TRUE)
    gene_expression_table(input$show_results[3])
  })
  output$relative_inter_4<-renderTable({
    source('global.R', local = TRUE)
    gene_expression_table(input$show_results[4])
  })
  output$relative_inter_5<-renderTable({
    source('global.R', local = TRUE)
    gene_expression_table(input$show_results[5])
  })
  output$relative_inter_6<-renderTable({
    source('global.R', local = TRUE)
    gene_expression_table(input$show_results[6])
  })
  output$relative_inter_1_plot<-renderPlot({
    source('global.R', local = TRUE)
    draw_plot(input$show_results[1],input$plottype)  
  })
  output$relative_inter_2_plot<-renderPlot({
   source('global.R', local = TRUE)
    draw_plot(input$show_results[2],input$plottype) 
  })
  output$relative_inter_3_plot<-renderPlot({
    source('global.R', local = TRUE)
    draw_plot(input$show_results[3],input$plottype) 
  })
  output$relative_inter_4_plot<-renderPlot({
    source('global.R', local = TRUE)
    draw_plot(input$show_results[4],input$plottype) 
  })
  output$relative_inter_5_plot<-renderPlot({
    source('global.R', local = TRUE)
    draw_plot(input$show_results[5],input$plottype) 
  })
  output$relative_inter_6_plot<-renderPlot({
    source('global.R', local = TRUE)
    draw_plot(input$show_results[6],input$plottype) 
  })
  output$plots<-renderUI({

  })
  
}

shinyApp(ui, server)
