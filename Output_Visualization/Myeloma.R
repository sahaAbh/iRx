##Chose the correct directory to work with. It should be "./Data_And_Codes/Output_Visualization
##----Specify input file details---
MatOut <- "MatOut_Myeloma.mat";
RdataOut <- "DataForFittingMyeloma.Rdata" ; 
##-----------------------------------
#-------Specify relevant directories------------
scriptsDir1 <- "../Scripts"
scriptsDir2 <- "../MetaData"
scriptsDir3 <-  "../Fig"
source(file.path(scriptsDir1, "PCR_ind_joint.R"))

##-- Load libraries---------------
##-- Install following packages before the running the code.
library(R.matlab); # For reading matlab files 
library(gplots);   # For barplots
library(ggplot2);
library(ggpubr); # needed for plotting p values
library(plyr)
library(scales)  # for using the function rescale in color palletting
library(reshape)
library(circlize)  # for circular mixed cluster representation 
library(pROC)
library(car) # for qqPlot
##----------------------------------

#-------------PCR methods------
rdata <- file.path(scriptsDir2, RdataOut)
pcr1 <- pcr_pred(rdata,nc=280,np=188)
pcr2 <- pcr_pred2(rdata,nc=280,np=188)

##--  Load the output data as obtained from MCMC procedure
Clst<- readMat(file.path(scriptsDir2, MatOut))
fp<-load(file.path(scriptsDir2, RdataOut));
Boxt<- list(IndResp= Clst$IndResp, predictp.dep= Clst$iRx.std, predictp.ind= Clst$Ni.std)
IndResp<- rep(0,188); 
IndResp[Boxt$IndResp==0]=1;IndResp[Boxt$IndResp==2]=2;
roc_obj <- roc(IndResp[IndResp!=2], Boxt$predictp.dep[IndResp!=2], direction="<")
auc(roc_obj)
roc_obj <- roc(IndResp[IndResp!=2], Boxt$predictp.ind[IndResp!=2], direction="<")
auc(roc_obj)
roc_obj <- roc(IndResp[IndResp!=2], pcr1$PCR_1[IndResp!=2], direction="<")
auc(roc_obj)
roc_obj <- roc(IndResp[IndResp!=2], pcr2$PCR_2[IndResp!=2], direction="<")
auc(roc_obj)
#   AUC:  PCR1 <- 39%, PCR2<- 41% , NI <- 46%  iRx <- 59%   
##----------------------------------
##-----------------------Figure 2(a) and 2(c): Violin plot and density plots 

#-------Figure 4(a)---[For final configuration of figure 2, use size: 11 vs 6]
#--------Labels for 4 methods
IndRespLevels1.1<- rep("Factor-Resp",188); IndRespLevels1.1[Boxt$IndResp==0]="Factor-Non-Resp";IndRespLevels1.1[Boxt$IndResp==2]="Factor-Others";
IndRespLevels1.2<- rep("Ind-Resp",188); IndRespLevels1.2[Boxt$IndResp==0]="Ind-Non-Resp";IndRespLevels1.2[Boxt$IndResp==2]="Factor-Others";
IndRespLevels1.3<- rep("PCR1-Resp",188); IndRespLevels1.3[Boxt$IndResp==0]="PCR1-Non-Resp";IndRespLevels1.3[Boxt$IndResp==2]="PCR1-Others";
IndRespLevels1.4<- rep("PCR2-Resp",188); IndRespLevels1.4[Boxt$IndResp==0]="PCR2-Non-Resp";IndRespLevels1.4[Boxt$IndResp==2]="PCR2-Others";
IndRespLevels1<- c(IndRespLevels1.1,IndRespLevels1.2,IndRespLevels1.3,IndRespLevels1.4)
#----------within method: labels for resp vs non-resp
IndRespLabels2.1<- rep("Resp",188);IndRespLabels2.1[Boxt$IndResp==0]="Non-Resp";IndRespLabels2.1[Boxt$IndResp==2]="Others";
IndRespLabels2<- rep(IndRespLabels2.1,4)

#----------over all 4 groups
Groups<- c(rep("Factor",188),rep("Independent",188), rep("PCR1",188), rep("PCR2",188));
plotgg1<-data.frame(Patients=IndRespLevels1,LogIC50=c(t(Boxt$predictp.dep),t(Boxt$predictp.ind), pcr1$PCR_1_std, pcr2$PCR_2_std ),Type=IndRespLabels2,Groups<-Groups)
plotgg2<-plotgg1[plotgg1$Type!="Others",]
my_comparisons <- list( c("Factor-Resp", "Factor-Non-Resp")  , c("Ind-Resp", "Ind-Non-Resp"), c("PCR1-Resp", "PCR1-Non-Resp"), c("PCR2-Resp", "PCR2-Non-Resp") )
set.seed(1234)
p3 <- ggplot(plotgg2,aes( Patients,LogIC50))
bold.10.text <- element_text(face = "bold", size = 10);
bold.text <- element_text(face = "bold");
p4<-p3+geom_violin(aes(fill = Type),show.legend=T,trim=F)+
  ylim(-4.5,6)+
  geom_jitter(height = 0, width = 0.1)+
  geom_boxplot(width=.1, outlier.colour=NA)+
  stat_compare_means(method = "t.test",method.args = list(alternative = "less"),comparisons = my_comparisons)+
  scale_x_discrete(labels = c("Non-resp","Resp","Non-resp","Resp","Non-resp","Resp","Non-resp","Resp"))+
  labs( title= "Myeloma Patients", x = "Patients' status", y = "Predicted LogIC50")+  
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title=element_text( hjust=0.5, vjust=0.5, face='bold'))+
#  scale_fill_manual(values= hue_pal()(2))+
  scale_fill_manual(values= grey.colors(5)[c(2,4)])+
  geom_text(x=1.5,y=4,label={"iRx"})+
  geom_text(x=3.5,y=4.5,label={"NI"})+
  geom_text(x=5.5,y=5.5,label={"PCR_IND"})+
  geom_text(x=7.5,y=6,label={"PCR_JOINT"})
#pdf(file=file.path(scriptsDir3,"Figure_2a_violin.pdf"), width=11, height=8.5);
pdf(file=file.path(scriptsDir3,"Figure_4a_violin_gray.pdf"), width=11, height=8.5);
#par(mfrow=c(2,1))
plot(p4)
dev.off()
#---- Figure 4(c): density plot indicating the predictive rule
#pdf(file=file.path(scriptsDir3,"Figure_2c_EP.pdf"));
pdf(file=file.path(scriptsDir3,"Figure_4c_EP_gray.pdf"));
d <- density(Boxt$predictp.dep) # returns the density data 
dtable <- data.frame(x=d$x ,y=d$y); 
ggplot(dtable, aes(x , y ))+
  geom_line(size = 1) +
#  geom_area(data = dtable[dtable$x < -0.5, ], fill = hue_pal()(2)[2])+
 # geom_area(data = dtable[dtable$x > 0.5, ], fill = hue_pal()(2)[1])+
  geom_area(data = dtable[dtable$x < -0.5, ], fill = gray.colors(5)[4])+
  geom_area(data = dtable[dtable$x > 0.5, ], fill = gray.colors(5)[2])+
  scale_x_continuous(name = "Standardized iRx scores",breaks = seq(-4, 4, 0.5))+
  theme(axis.text.x = bold.10.text,plot.title=element_text( hjust=0.5, vjust=0.5, face='bold'))+
  geom_text(x=3,y=0.2,label={"EP=0.99"})+
  labs( title= "Myeloma Patients", y = "Density")
##------------------------------------------------
dev.off();

##---------------------------Figure 3(a) and 3(b)
#----------------------------- computing more cluterspecific details------------
#--- Collecting the names of the samples, cluster indicator and sensitivity status and sensitivity values
idx4<- Clst$idx
idx4<- as.data.frame(idx4);
#--v4 indicates sample no.
#--v2 1 for cell line and 0 for patients
#--v1 cluster no.
#--v3 indicator for cell line(3) or Patients responder(1) or Non responder(0) or Others(2)
CellCnt<-aggregate(V2~V1,data=idx4 ,"sum") # Gives total number of cell lines
TotCnt<-aggregate(V2~V1,data=idx4 ,"length") # Gives total size
#--Without the patient sub-classification
idx5Sum<-cbind(TotCnt,V3=CellCnt[,2],V4=TotCnt[,2]-CellCnt[,2]); #---v4 no of patient
idx5Sum<-idx5Sum[order(idx5Sum$V2,decreasing = T),]
idx4.1<- idx4[idx4[,3]!=2,] # patients removed who are neither responder nor non-responder
f1<-tapply(idx4.1$V3,idx4.1$V1, function(x)sum(x==3))
f2<-tapply(idx4.1$V3,idx4.1$V1, function(x)sum(x==1))
f3<-tapply(idx4.1$V3,idx4.1$V1, function(x)sum(x==0))
TotCnt<-aggregate(V2~V1,data=idx4.1 ,"length") # Gives total size
idx5Sum.1<-cbind(TotCnt,V3=f1,V4=f2,V5=f3);
idx5Sum.1<-idx5Sum.1[order(idx5Sum.1$V2,decreasing = T),]
IndRespLevels<- rep("Responder",188)  
IndRespLevels[Boxt$IndResp==0]="Non-Responder";
IndRespLevels[Boxt$IndResp==2]="Others";
plotgg<-data.frame(Patients=IndRespLevels,Dep=t(Boxt$predictp.dep),ind=t(Boxt$predictp.ind))
idx4.2 <- idx4[order(idx4[,4]),]
#rownames(idx4.2)<- Allnames_CP;
#idx4.2 <- idx4[order(idx4[,4]),]
#rownames(idx4.2)<- Allnames_CP;
#-------------------------------------------------
#----------------------------UNCOMMENT THE CODE BELOW TO LIST THE CELL LINES  
NameId<-  Allnames_CP;
NameId[idx4.2$V2==0]= "";
# UNCOMMENt BELOW TO STORE THE AVATARS
#tmp<- idx4.2[idx4.2$V2==1, c("V1","V4")];
#tmp$V4 <- NameId[tmp$V4];
idx5.2 <- tapply(NameId,idx4.2$V1,function(x){paste(x, sep="",collapse=",")})
TotCnt<-aggregate(V2~V1,data=idx4.1 ,"length")
idx5Sum.2<-cbind(TotCnt,V3=idx5.2);
idx5Sum.2<- idx5Sum.2[order(idx5Sum.2$V2,decreasing = T),] 
#Map<- cbind(old=idx5Sum.2[,1],new=1:11); Map=Map[order(Map[,1]),];
#tmp$Cl= Map[tmp$V1,"new"];
#write.csv(tmp, file="ClusterAvatarsBort.csv", row.names=F)
idx5Sum.2$V1<- sort(idx5Sum.2$V1); 
colnames(idx5Sum.2)=c("Cluster index","Cluster size", "Cell lines");
write.csv(idx5Sum.2, file="ClusterCelllinesBort.csv", row.names=F);
#-------------------------------------------------------------------------
idx4_sensitivity <- cbind(idx4.2, sensitivity=0);
idx4_sensitivity[idx4_sensitivity$V3!=3,"sensitivity"]= plotgg$Dep;
#------------choosing mixed cluster------------------------------------
#---idX6$V2 V1 has the mixed cluster information 
#----Mixed clusters 
idx5_mixed<- idx5Sum.1[as.logical((idx5Sum.1$V3 >0) * ((idx5Sum.1$V2-idx5Sum.1$V3) >0)), ];
#---------------------Only iRx score------------------------------------
Cl=idx5Sum.1[,1]; #Mixed <- paste0('G',Cl[c(1:3,5:6,9:10)]); 
Allclusters <- paste0('G',Cl); #Cl1= c(1:3,5:6,9:10); 
idx4_patients_all= idx4_sensitivity[idx4_sensitivity$V3!=3,];
idx4.1 = idx4_sensitivity[idx4_sensitivity$V3!=2,]
idx4_patients= idx4.1[idx4.1$V3!=3,];
scale <- sd(idx4_patients_all$sensitivity)/2
Mn <- mean(idx4_patients_all$sensitivity)
Cutoff_U=scale+Mn;
Cutoff_L=-1*scale+Mn;

S1= 1*(idx4_patients$sensitivity > Cutoff_U); 
S1[idx4_patients$sensitivity < Cutoff_L ]=-1;
#----------------------------
Data1=cbind(idx4.1 , Sen= 2);
Data1[Data1$V3!=3,"Sen"]=S1;
TotCnt<-aggregate(V2~Sen+V1,data=Data1,"length") # Gives total size
Data2<-reshape(as.data.frame(TotCnt), direction="wide", v.names = 'V2', timevar = "Sen",
               idvar = "V1")
rownames(Data2)= paste0('G',Data2[,1]);
Data3= Data2[Allclusters,];
Data3[is.na(Data3)]=0;
names.arg= paste0('C',1:11);
pdf(file=file.path(scriptsDir3,"Figure_5_chord.pdf"), width=16, height=8.5);
par(mfrow=c(1,2))
barplot2(t(as.matrix(Data3[,-1])), beside = F ,horiz=F, names.arg=names.arg,xlab="Cluster Index",main= "Split by Cluster", ylab="Frequency",col=c("blue","yellow","beige","red"))
legend("topright",legend=c("Responder","Neutral","Non-responder","Cell lines"),fill=c("blue","yellow","red","beige"))
#---------------------------3(b)-----------------
#========================================  chord diagram===========================
idx4.1[ idx4.1$V3!=3,"V3"]=S1;
idx4.1[idx4.1$V3==3,"V3"]=2;
RespTag=c("Responder", "Neutral", "Non-responder")
# requires circlize library:  
#----idx4.1 MUST have V1 (CLUSTER INDEX), V2(PATIENT INDEX), V3(RESPONDER INDEX), V4(SAMPLE INDEX)
#----idx5Sum.1  MUST have V1 (CLUSTER INDEX), V2(TOTAL SAMPLE SIZE), V3(CELL LINE SIZE), V4(PATIENT RESP SIZE), V4(PATIENT NON-RESP SIZE) 
PCONLY = idx5Sum.1[ as.logical((!((idx5Sum.1$V4==0) * (idx5Sum.1$V5==0))) * (idx5Sum.1$V3 >0)), 1 ]  
idx4.PC = idx4.1[idx4.1$V1 %in% PCONLY ,];
idx4.C = idx4.PC[idx4.PC$V2==1,c("V1","V4")]; idx4.C<-idx4.C[order(idx4.C$V1,idx4.C$V4),]
idx4.P = idx4.PC[idx4.PC$V2==0,c("V1","V3","V4")]; idx4.P<-idx4.P[order(idx4.P$V1,idx4.P$V4),]
Adj<-merge(idx4.C,idx4.P,by="V1",suffixes=c(".C",".P"))
Rnk<- cbind(Rank=1:length(PCONLY),V1=PCONLY)
Adj1<-merge(Rnk,Adj,by="V1");Adj1<-Adj1[order(Adj1$Rank,Adj1$V4.C,Adj1$V4.P),]
#df = data.frame(from = AllNames[Adj1$V4.C],
#                to = RespTag[Adj1$V3+1],
#                value = 1,
#                stringsAsFactors = FALSE)
Adj2<-aggregate(V4.P~V1+Rank+V4.C+V3,data=Adj1,length);Adj2<-Adj2[order(Adj2$Rank,Adj2$V4.C,Adj2$V3),]
df1 = data.frame(from = Allnames_CP[Adj2$V4.C],
                 to = RespTag[Adj2$V3+2],
                 value = Adj2$V4.P,
                 stringsAsFactors = FALSE)
ClSize<-idx5Sum.1[as.character(PCONLY),"V3"]
circos.par(start.degree = 180,gap.after = c(rep(0,ClSize[1]-1),5,rep(0,ClSize[2]-1),5,rep(0,ClSize[3]-1),5,rep(0,ClSize[4]-1),5,
                                            rep(0,ClSize[5]-1),5,rep(0,ClSize[6]-1),5,rep(0,ClSize[7]-1),5,rep(0,ClSize[8]-1),5,5,5,5),clock.wise = FALSE)
chordDiagram(df1, grid.col=c(rep("grey",sum(ClSize)),c("blue","yellow","red")),col = Adj2$Rank,link.border= Adj2$Rank, annotationTrack = c("grid",axis),transparency=0.2, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df1))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if((xplot[2] - xplot[1]) < -20) {
    circos.text(mean(xlim), mean(ylim), sector.name, facing = "outside",
                niceFacing = TRUE,  col = "black")
    print((xplot[2] - xplot[1]))
  }else{
    if((xplot[2] - xplot[1]) < -2) {
      circos.text(mean(xlim), mean(ylim), sector.name, facing = "clockwise",
                  niceFacing = TRUE,  col = "black")
      print((xplot[2] - xplot[1]))
    }
  } 
}, bg.border = NA)
circos.clear()
##-------------------------------------------------------------------------
dev.off();
##-----------------------------------------------------------------------------
LambdaMat<-Clst$L.mh
GenFilter<-function(Col,Percentile){FctOrder= order(LambdaMat[,Col],decreasing = T);Fct<-LambdaMat[FctOrder,Col];Lim<-which(cumsum(Fct)>=Percentile)[1];Topgenes=FctOrder[1:Lim];}
Genes<- GenFilter(1,0.95)
for(i in 2:11)
  {
  Genes <- union(Genes,GenFilter(i,0.95))
  print(length(Genes))
  }
#-----------Total number of top genes 978 

LambdaMat1<- LambdaMat[Genes,]
Mycol<-colorRampPalette(c("black","blue", "green","yellow","red"),space="Lab")(5)
#------------Breaks in lojg10 scale-------------
col_breaks =c(min(log10(LambdaMat1))-0.5, quantile(log10(LambdaMat1),c(0.25,0.5,0.75,0.95)),max(log10(LambdaMat1))+0.5)
pdf(file = NULL);
hm<-heatmap.2(log10(LambdaMat1),
              Rowv=T,
              #cellnote = LambdaMat1,same data set for cell labels
              main = "Gene clusters across factors", # heat map title
              # notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              #margins =c(1,1),     # widens margins around plot
              breaks=col_breaks,
              xlab= "Factors",
              ylab= "Genes",
              Colv=F,
              col=Mycol ,
              dendrogram="none",
              scale= "none",
              labRow="",
              labCol= paste("Fact-",1:11,sep=""),
              key.xlab="log10(Value)");
dev.off();
LambdaMat2<- as.data.frame(log10(LambdaMat1[rev(hm$rowInd),]))
LambdaMat2$Name1 <- as.factor(1:nrow(LambdaMat2))
LambdaMat2.m<- melt(LambdaMat2, id.vars= "Name1")
LambdaMat2.m$var2 = as.numeric(LambdaMat2.m$variable) + 5

# Labels and breaks need to be added with scale_y_discrete.
y_labels =  paste("F",c(1:11),sep="")
y_breaks = c(6:16)

library(scales)
p2 = ggplot(LambdaMat2.m, aes(x=Name1, y=var2, fill=value)) +
  geom_tile()+  
  #scale_fill_gradient(low = "white", high = "steelblue",breaks=col_breaks3) +
  scale_fill_gradientn(colours=Mycol, values=rescale(col_breaks),
  name="Value(log10)",breaks=c(-21.65,-13.85,-5.00,-2.31),
  labels=c("-21.65","-13.85","-5.00","-2.31")) +
  # original scale: c("2.23e-22","1.41e-14","1e-05","4.8e-03")
  #ylim(c(0, 20)) +
  scale_y_continuous(breaks=y_breaks,labels= y_labels,limits=c(0,17))+
  coord_polar(theta="x") +
  theme(panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=5))
pdf(file=file.path(scriptsDir3,"Figure_6a_loadings.pdf"), width=11, height=8.5);
plot(p2)
dev.off()
#-----------------------------------------------------------------------------
fp<-load(file.path(scriptsDir2,"DataForFittingMyeloma.Rdata"))
pdf("../Fig/myelomaqq.pdf",width=5,height=5)
qqPlot(D, id=F, main="Myeloma: QQ plot", ylab= "transformed IC50", col.lines="black")
dev.off()
