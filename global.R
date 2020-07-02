# functions used in app.R
# read the input file
read_data<-function(input){
  file2table<-read_excel(input, col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[-(1:37),]
  colnames(file2table)<-read_excel(input, col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[37,]
  file2table<-data.frame(file2table)[,c(1,3,4,8,9,10,15,23)]
  return(file2table)
}
# get all sample data
get_all_sample_data <- function (readfile){
  all_sample_data<-NULL
  for (sample_name in row.names(table(readfile[,'Sample.Name']))){
    sample_A<-readfile[readfile$Sample.Name==sample_name,]
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
}

# get relative samples data
get_relative_samples_data<-function(all_sample_data){
  relative_samples_data<-NULL
  for (sample_name in row.names(table(all_sample_data$Sample.Name))){
    A_sample_data<-all_sample_data[all_sample_data$Sample.Name==sample_name,]
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
}

get_Ctrls.Mean.dCt_table <- function (relative_data_intra){
  Samples.Name<-row.names(table(relative_data_intra$Sample.Name))
  genes<-row.names(table(relative_data_intra$Target.Name))
  Ctrl.Samples<-input$controls
  Ctrls.Mean.dCt<-NULL
  for (gene in genes ) {
    Ctrls.dCt<-NULL
    for (Ctrl.Sample in Ctrl.Samples){
      Ctrl.dCt<-relative_data_intra[relative_data_intra$Sample.Name==Ctrl.Sample,][relative_data_intra[relative_data_intra$Sample.Name==Ctrl.Sample,]$Target.Name==gene,]$dCt
      Ctrls.dCt<-c(Ctrls.dCt,Ctrl.dCt)
    }
    Ctrl.Mean.dCt<-mean(Ctrls.dCt)
    Ctrls.Mean.dCt<-c(Ctrls.Mean.dCt,Ctrl.Mean.dCt)
  }
  Ctrls.Mean.dCt_table<-data.frame(genes)
  Ctrls.Mean.dCt_table<-cbind(Ctrls.Mean.dCt_table,Ctrls.Mean.dCt)
  colnames(Ctrls.Mean.dCt_table)<-c('Target.Name','dCt')
  return(Ctrls.Mean.dCt_table)
}
# set gene and treatment groups
get_relative_data_inter<-function(relative_data_intra){
  relative_data_inter<-relative_data_intra
  ddCt_list<-NULL
  # sorting gene group
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
  # sorting treatment group
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
}

gene_expression_table <- function (input){
  if (is.empty(input)){
    return(NULL)
  } else
    return(relative_data_inter()[relative_data_inter()$Target.Name==input,])
}

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
# draw plots
box_plot<-function(input){
  p<-ggplot(data = relative_data_inter()[relative_data_inter()$Target.Name==input,], aes(x=Group, y=exprs)) +
    geom_boxplot() +
    ggtitle(label = input) +
    scale_x_discrete(name ="Groups", 
                     limits= rev(row.names(table(relative_data_inter()[relative_data_inter()$Target.Name==input,]$Group))) )+
    stat_summary(fun.data=data_summary, 
                 geom='pointrange', color="red")
  return(p)
}

dot_plot<-function(input){
  p<-ggplot(data = relative_data_inter()[relative_data_inter()$Target.Name==input,], aes(x=Group, y=exprs)) +
    geom_count() +
    ggtitle(label = input) +
    scale_x_discrete(name ="Groups", 
                     limits= rev(row.names(table(relative_data_inter()[relative_data_inter()$Target.Name==input,]$Group))) )+
    stat_summary(fun.data=data_summary, 
                 geom='pointrange', color="red")
  return(ggplotly(p))
}

dot_box_plot<-function(input){
  p<-ggplot(data = relative_data_inter()[relative_data_inter()$Target.Name==input,], aes(x=Group, y=exprs)) +
    geom_boxplot() +
    geom_count()  +
    ggtitle(label = input) +
    scale_x_discrete(name ="Groups", 
                     limits= rev(row.names(table(relative_data_inter()[relative_data_inter()$Target.Name==input,]$Group))) )+
    stat_summary(fun.data=data_summary, 
                 geom='pointrange', color="red")
  return(p)
}

draw_plot <- function(input,plottype) {
  if (plottype=='Dot-plot'){
    dot_plot(input)
  } else if (plottype=='Box-plot'){
    box_plot(input)
  } else if (plottype=='Dot-Box-plot'){
    dot_box_plot(input)
  }
}
