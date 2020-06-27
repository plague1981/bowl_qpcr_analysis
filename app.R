library(shiny)
library(shinythemes)
library(shinydashboard)
library(readxl)
library(xlsx)
library(ggplot2)
library(rapportools)

ui<- dashboardPage(
     dashboardHeader(title = 'Bowl\'s project'),
     dashboardSidebar(
       fileInput(inputId = 'file',label = 'Select input file:',multiple = FALSE),
       tableOutput('filename'),
       actionButton(inputId = "analyze", label = "Analyze")
     ),
     dashboardBody(
       navbarPage(title ='Data analysis', 
         tabPanel('Raw data', icon = icon('file'),
           tableOutput('rawdata')
         ),
         tabPanel('Calibration', icon = icon('file'),
           tableOutput('calibration')
         )
             
       )
           
     )
)


server <- function(input, output, session){
  # functions
  readfile<- eventReactive(input$analyze,{ 
    file2table<-read_excel(input$file$datapath, col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[-(1:37),]
    colnames(file2table)<-read_excel(input$file$datapath, col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[37,]
    file2table<-data.frame(file2table)[,c(1,3,4,8,9,10,15,23)]
    return(file2table)
  })
  calibration<- eventReactive(input$analyze,{


  })
  output$filename<-renderTable({
    input$file[1]
  })
  output$rawdata<-renderTable({
    readfile()
  })
  output$calibration<-renderTable({
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
}

shinyApp(ui, server)