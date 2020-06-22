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


sample_A<-newdata[newdata$Sample.Name=='327 Undiffer',]
hq_data<-NULL
for (gene in rownames(table(sample_A$Target.Name))){
  if (sample_A[sample_A$Target.Name==gene, ]$Ct.Threshold[1]<=0.1){
    hq_data<-rbind(hq_data, sample_A[sample_A$Target.Name==gene, ][1,c(2,3,5)])
  } 
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
    if (which.min(c(a,b,c))==a){
      new_Ct.Mean<-(as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[1])+as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[2]))/2
      new_Ct.Mean_list<-c(new_Ct.Mean_list,new_Ct.Mean)
    } else if (which.min(c(a,b,c))==b){
      new_Ct.Mean<-(as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[2])+as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[3]))/2
      new_Ct.Mean_list<-c(new_Ct.Mean_list,new_Ct.Mean)
    } else {
      new_Ct.Mean<-(as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[1])+as.numeric(sample_A[sample_A$Target.Name==gene,]$CT[3]))/2
      new_Ct.Mean_list<-c(new_Ct.Mean_list,new_Ct.Mean)
    }
    lq_gene<-c(lq_gene,gene)
  } 
}
