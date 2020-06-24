library(readxl)

dirPath<-'/home/joey/Documents/WanYing'
setwd(dirPath)

tryCatch(
  expr = {
    dir.create(file.path(getwd(), 'Results'))
    message("Successfully created 'results' folder." )
  },
  error = function(e){
    message('Caught an error!')
    print(e)
  },
  warning = function(w){
    message('Caught an warning!')
    print(w)
  },
  finally = {
    message('All set')
  }
)
# Set export directory
export_dir<-paste0(dirPath,'/Results')

files<-list.files(path = getwd(), pattern = 'xls')

rawdata<-read_excel(files[1], col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[-(1:37),]
col_names<-read_excel(files[1], col_types="text", range=cell_cols(c("B","X")), col_names = TRUE ,trim_ws = TRUE, sheet = 'Results')[37,]
colnames(rawdata)<-col_names
newdata<-data.frame(rawdata)[,c(1,3,4,8,9,10,15,23)]

all_sample_data<-NULL
for (sample_name in row.names(table(newdata[,'Sample.Name']))){
  sample_A<-newdata[newdata$Sample.Name==sample_name,]
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

# differential_gene_express in samples

relative_samples_data<-NULL
for (sample_name in row.names(table(all_sample_data$Sample.Name))){
  A_sample_data<-all_sample_data[all_sample_data$Sample.Name==sample_name,]
  Crtl_Ct<-A_sample_data[A_sample_data$Target.Name=='GAPDH',]$Ct.Mean
  genes<-NULL
  genes_dCt<-NULL
  relative_sample_data<-NULL
  for (gene in A_sample_data$Target.Name){
    if (gene!='GAPDH'){
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

# list sample names
Samples.Name<-row.names(table(relative_samples_data$Sample.Name))
Ctrl.Samples<-c('wt_272 Undiffer','wt_279 Undiffer','wt_328 Undiffer','wt_339 Undiffer')
Ctrls.Mean.dCt<-NULL
for (gene in genes ) {
  Ctrls.dCt<-NULL
  for (Ctrl.Sample in Ctrl.Samples){
    Ctrl.dCt<-relative_samples_data[relative_samples_data$Sample.Name==Ctrl.Sample,][relative_samples_data[relative_samples_data$Sample.Name==Ctrl.Sample,]$Target.Name==gene,]$dCt
    Ctrls.dCt<-c(Ctrls.dCt,Ctrl.dCt)
  }
  Ctrl.Mean.dCt<-mean(Ctrls.dCt)
  Ctrls.Mean.dCt<-c(Ctrls.Mean.dCt,Ctrl.Mean.dCt)
}
Ctrls.Mean.dCt_table<-data.frame(genes)
Ctrls.Mean.dCt_table<-cbind(Ctrls.Mean.dCt_table,Ctrls.Mean.dCt)
colnames(Ctrls.Mean.dCt_table)<-c('Target.Name','dCt')

# Create ddCt table


ddCt_list<-NULL
for (n in 1:nrow(relative_samples_data)){
  if (relative_samples_data[n,]$Target.Name=='BSP'){
    ddCt<-(relative_samples_data[n,'dCt']-Ctrls.Mean.dCt_table[Ctrls.Mean.dCt_table$Target.Name=='BSP',]$dCt)
  } else if (relative_samples_data[n,]$Target.Name=='Col1a1') {
    ddCt<-(relative_samples_data[n,'dCt']-Ctrls.Mean.dCt_table[Ctrls.Mean.dCt_table$Target.Name=='Col1a1',]$dCt)
  } else if (relative_samples_data[n,]$Target.Name=='Col2a1') {
    ddCt<-(relative_samples_data[n,'dCt']-Ctrls.Mean.dCt_table[Ctrls.Mean.dCt_table$Target.Name=='Col2a1',]$dCt)
  } else if (relative_samples_data[n,]$Target.Name=='OSX') {
    ddCt<-(relative_samples_data[n,'dCt']-Ctrls.Mean.dCt_table[Ctrls.Mean.dCt_table$Target.Name=='OSX',]$dCt)
  } else if (relative_samples_data[n,]$Target.Name=='Runx2') {
    ddCt<-(relative_samples_data[n,'dCt']-Ctrls.Mean.dCt_table[Ctrls.Mean.dCt_table$Target.Name=='Runx2',]$dCt)
  } else if (relative_samples_data[n,]$Target.Name=='SOX9') {
    ddCt<-(relative_samples_data[n,'dCt']-Ctrls.Mean.dCt_table[Ctrls.Mean.dCt_table$Target.Name=='SOX9',]$dCt)
  }
  ddCt_list<-c(ddCt_list,ddCt)
}

exprs<-2^-ddCt_list
relative_samples_data$ddCt<-ddCt_list
relative_samples_data$exprs<-exprs
wt_undiff<-NULL

wt_undiffs<-wt_diffs<-mu_undiffs<-mu_diffs<-NULL
for (n in 1:nrow(relative_samples_data[relative_samples_data$Target.Name=='BSP',])){
  if (grepl('^wt', relative_samples_data[relative_samples_data$Target.Name=='BSP',]$Sample.Name[n]) & grepl('Undiffer$', relative_samples_data[relative_samples_data$Target.Name=='BSP',]$Sample.Name[n])){
    wt_undiff<-relative_samples_data[relative_samples_data$Target.Name=='BSP',][n,]
  }
  wt_undiffs<-rbind(wt_undiffs,wt_undiff)
}

