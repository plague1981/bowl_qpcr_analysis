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
    geom_dotplot(binaxis='y', stackdir='center',) +
    ggtitle(label = input) +
    scale_x_discrete(name ="Groups", 
                     limits= rev(row.names(table(relative_data_inter()[relative_data_inter()$Target.Name==input,]$Group))) )+
    stat_summary(fun.data=data_summary, 
                 geom='pointrange', color="red")
  return(p)
}

dot_box_plot<-function(input){
  p<-ggplot(data = relative_data_inter()[relative_data_inter()$Target.Name==input,], aes(x=Group, y=exprs)) +
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center',) +
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



