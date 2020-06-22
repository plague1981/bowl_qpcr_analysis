dirPath<-'C:/Users/Changyi.Lin/Documents/bowl'
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

A_sample_data<-all_sample_data[all_sample_data$Sample.Name=='mu_270 Differ',]
Crtl_Ct<-A_sample_data[A_sample_data$Target.Name=='GAPDH',]$Ct.Mean
genes<-NULL
genes_de<-NULL
for (gene in A_sample_data$Target.Name){
  if (gene!='GAPDH'){
  Gene_Ct<-A_sample_data[A_sample_data$Target.Name==gene,]$Ct.Mean
  gene_de<-2^(-(as.numeric(Gene_Ct)-as.numeric(Crtl_Ct)))
  genes<-c(genes,gene)
  genes_de<-c(genes_de,gene_de)
  }
}
barplot(genes_de, names.arg = genes, cex.names=0.8)
