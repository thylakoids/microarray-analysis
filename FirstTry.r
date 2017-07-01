#process zhu's data 
#first try
rm(list=ls())
#wap function to handle constant 
t.test.constant=function(x,y){
  tryCatch( return(t.test(x,y,var.equal = T)$p.value),
            error = 
              function(e){return(1)})}
#read data(already normalizd)
data.all=read.csv('startingData.csv',row.names = 1)
#three part
data.raw=data.all[,1:10]
data.log2normalized=data.all[,11:20]
data.flags=data.all[,21:30]
#p-value
data.pvalue=apply(data.log2normalized,1,function(x)
  {t.test.constant(x[c(6,8)],x[c(4,5,10)])})
#ap
AP=apply(data.flags[c(4,5,6,8,10)],1,paste,collapse='')
genes.present = names(AP[AP != "AAAAA"])
#mean
data.sensitivity.mean=apply(data.log2normalized[c(6,8)],1,mean)
data.resistence.mean=apply(data.log2normalized[c(10,4,5)],1,mean)
data.sensitivity.mean.present=data.sensitivity.mean[genes.present]
data.resistence.mean.present=data.resistence.mean[genes.present]
#ratio
data.ratio=apply(data.log2normalized,1,function(x)
  {log2y=mean(x[c(6,8)]-mean(x[c(10,4,5)]));2**log2y})
data.ratio.present=data.ratio[genes.present]
data.regulation=sapply(data.ratio,function(x){ifelse(x>1,'up','down')})
data.regulation.present=data.regulation[genes.present]
#fdr
data.present=data.log2normalized[genes.present,c(6,8,10,4,5)]
data.pvalue.present=data.pvalue[genes.present]
data.fdr.present=p.adjust(data.pvalue.present,method = 'fdr')
#combine
data.plus=cbind(data.present,data.sensitivity.mean.present,
                data.resistence.mean.present,
                data.pvalue.present,
                data.fdr.present,data.ratio.present,
                data.regulation.present)
#DE
genes.DE=genes.present[data.pvalue.present<0.01]
data.DE.plus=data.plus[genes.DE,]




#volcano plot
library(ggplot2)
library(RColorBrewer)
Set3=brewer.pal(3,'Set3')
labels=c('not differential expresssed','up regulated','down regulated')
threshold1=as.numeric(data.ratio.present>2&data.pvalue.present<0.05)
threshold2=as.numeric(data.ratio.present<0.5&data.pvalue.present<0.05)
data.plus$threshold=labels[(threshold1+2*threshold2)+1];

g=ggplot(data=data.plus,
         aes(x=log(data.ratio.present,2),
         y=-log(data.pvalue.present,10),
         color=threshold))+
  scale_color_manual(values =c('green','black','red'))+
  geom_point(alpha=0.4,size=3)+
  geom_vline(xintercept = c(-1,1),color='green',size=1.3)+
  geom_hline(yintercept = -log(0.05,10),color='green',size=1.3)+
  xlim(c(-10,10))+ylim(c(0,8.3))+
  xlab('log2 fold change')+ylab('-log10 p-value')+
  theme(legend.position = "top")+
  coord_fixed(ratio=2) 
         
g

#boxplot
library(reshape2)
library(ggplot2)
colnames(data.log2normalized)=paste('ESCC.',1:10,sep='')
df.data.log2normalized=melt(data.log2normalized)
g=ggplot(df.data.log2normalized,
       aes(x=variable,
           y=value)
       )+
  geom_boxplot(fill=rgb(0,1,1),color=rgb(0,0,1))+
  xlab('Sample')+ylab('Normalized Intensity Values')+
  scale_y_continuous(breaks=seq(2,20,2),minor_breaks = NULL)+
  theme(axis.text.x=element_text(angle=45,vjust = 0.2,hjust=0.5))+
  coord_fixed(ratio=1/2) 
g

#scatter plot
threshold1=as.numeric(data.ratio.present>2)
threshold2=as.numeric(data.ratio.present<0.5)
data.plus$threshold.scatter=labels[(threshold1+2*threshold2)+1];
g=ggplot(data.plus,
         aes(y=data.sensitivity.mean.present,
             x=data.resistence.mean.present,
             color=threshold.scatter
             )
         )+
  scale_color_manual(values =c('green','black','red'))+
  geom_point(alpha=0.3,size=2.5)+
  geom_abline(slope=1,intercept = 1,size=1.3,color='gray',
              linetype='dashed')+
  geom_abline(slope=1,intercept =(-1),size=1.3,color='gray',
              linetype='dashed')+
  coord_fixed(ratio=1)
g
#heat map
library(gplots)
heatmap.2(as.matrix(data.present)[genes.DE,],
        col=greenred(100),
        ColSideColors = c(rep('red',2),rep('blue',3)),
        cexRow = 0.3,
        scale = 'none',
  #      density.info = 'none',
        symkey=FALSE,
        trace="none"
        )
