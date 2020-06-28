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

dot_plot<-function(input){
  p<-ggplot(data = relative_data_inter()[relative_data_inter()$Target.Name==input,], aes(x=Group, y=exprs)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  ggtitle(label = input) +
  scale_x_discrete(name ="Groups", 
                   limits= rev(row.names(table(relative_data_inter()[relative_data_inter()$Target.Name==input,]$Group))) )+
  stat_summary(fun.data=data_summary, 
               geom="pointrange", color="red")
  return(p)
}
