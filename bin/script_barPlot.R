
# 1. The dinucleotide plot. The proportion of each dinucleotide is weighted by the sum of expression score of dominant TSSs that with such dinucleotide. 
### I DID NOT subset the TCs mapping to promoters in this case, but it's achievable since we generated tmp as GRanges before plotting, thus we could always subset the dataframe with TCs annotated with promoter at the start.
library(viridis)
bsg<-BSgenome.Hsapiens.UCSC.hg38 #input related Bs.genome obj

for (i in 1:length(sampleLabels(CAGEr_object))) {
  print(i)
  tmp<-as.data.frame(tagClustersGR(CAGEr_object,
                                   sample = sampleLabels(CAGEr_object)[i],
                                   qLow = 0.1,qUp = 0.9))
  
  
  tmp<-GRanges(seqnames = tmp$seqnames,
                ranges =  
                  IRanges(start = tmp$dominant_ctss.pos,
                          end = tmp$dominant_ctss.pos),
                strand = tmp$strand,
                score=tmp$dominant_ctss.score,
                seqlengths = seqlengths(bsg))
  
  
  tmp<-promoters(tmp,
                 upstream = 1,downstream = 1)
  
  tmp<-tmp[width(trim(tmp)) == 2]
  
  tmp$dinucleotide<-as.data.frame(getSeq(bsg,tmp))$x
  
  A<-as.data.frame(tmp)[,c("score","dinucleotide")]
  dinucleotide<-unique(A$dinucleotide)
  B<-data.frame(dinucleotide=dinucleotide,sum_score=NA)
  for (l in dinucleotide) {
    
    B[which(B$dinucleotide==l),]$sum_score<-sum(A[which(A$dinucleotide==l),]$score)
  
  }
  
  B$proportion<-B$sum_score/sum(B$sum_score)
  
  B<-B[order(B$proportion,decreasing=F),]
  B$dinucleotide<-factor(B$dinucleotide,levels = B$dinucleotide)
  
  ggplot(data=B,aes(x=factor(dinucleotide),
                       y=(proportion*100),fill=dinucleotide))+
  geom_bar(stat="identity",
           colour = "black", 
           size = 0.3,
           position=position_dodge())+
  ggtitle(paste0(sampleLabels(CAGEr_object)[i], ": dominant TSS, -/+ 1bp position,
                 proportion weighted by the sum of dominant TSS score"))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=20),
        legend.key.size = unit(0.8,"cm"))+
  labs(fill="")+
  scale_fill_viridis(option = "magma",discrete = TRUE)
  
 ggsave(filename = paste0(sampleLabels(CAGEr_object)[i],"_DTSSs_weighted.pdf"),
        width = 12,
        height=10)

}



### 2. Alternatively, there is a script for generating dinucleotide barplot based on the Quantiles of expression of dominant peaks of Tag clusters, which might be considered to be included in the pipeline:
## What it does: plot dinucleotide distribution around dominant TSSs of TCs in each sample, individually based on dominant TSS score: top quantile 100%-75%, middle quantile 75%-25%, bottom quantile 25%-0%
library(viridis)

 #### For-loop for QC of each samples
for (i in 1:length(sampleLabels(CAGEr_object))) {
  print(i)
  
  
  tmp<-as.data.frame(tagClustersGR(CAGEr_object,
                                   sample = sampleLabels(CAGEr_object)[i],
                                   qLow = 0.1,qUp = 0.9))
  quantile25<-as.numeric(quantile(tmp$dominant_ctss.score))[2]
  quantile50<-as.numeric(quantile(tmp$dominant_ctss.score))[3]
  quantile75<-as.numeric(quantile(tmp$dominant_ctss.score))[4]
  
  quantileValue<-c(quantile25,quantile50,quantile75)
  
  tmp2<-GRanges(seqnames = tmp[which(tmp$dominant_ctss.score<quantileValue[1]),]$seqnames,
                ranges =  
                  IRanges(start = tmp[which(tmp$dominant_ctss.score<quantileValue[1]),]$dominant_ctss.pos,
                          end = tmp[which(tmp$dominant_ctss.score<quantileValue[1]),]$dominant_ctss.pos),
                strand = tmp[which(tmp$dominant_ctss.score<quantileValue[1]),]$strand,
                seqlengths = seqlengths(bsg))
  tmp2<-promoters(tmp2,upstream = 1,downstream = 1)
  tmp2<-tmp2[width(trim(tmp2)) == 2]
  
  c<-melt(table(getSeq(bsg,tmp2)))
 # c<-c[c$Var1%notin%c("NN","NG"),]
  c<-c[!is.na(c$Var1),]
  
  d<-mutate(c,perc=value/sum(value))
  d$value<-NULL
  
  d<-d[order(d$perc,decreasing=F),]
  d$Var1<-factor(d$Var1,levels = d$Var1)
  
  ggplot(data=d,aes(x=factor(Var1),
                    y=(perc*100),
                    fill=Var1))+
    geom_bar(stat="identity",
             colour = "black", 
             size = 0.3,
             position=position_dodge())+
    ggtitle(paste0("summary of Dominant_TSS_1st Quant, ",
                   sampleLabels(CAGEr_object)[i]," , -/+ 1 position"))+
    coord_flip()+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size=20),
          legend.key.size = unit(0.8,"cm"))+
    labs(fill="")+
    scale_fill_viridis(option = "magma", discrete = TRUE)
  
  ggsave(filename = paste0("./Dominant_TSS_1stQuant_",sampleLabels(CAGEr_object)[i],".pdf"),width = 12,height=10)
  
  tmp3<-GRanges(seqnames = tmp[which(tmp$dominant_ctss.score>quantileValue[1] & tmp$dominant_ctss.score<quantileValue[3] ),
  ]$seqnames,
  ranges = 
    IRanges(start = tmp[which(tmp$dominant_ctss.score>quantileValue[1] & tmp$dominant_ctss.score<quantileValue[3] ),
    ]$dominant_ctss.pos,
    end = tmp[which(tmp$dominant_ctss.score>quantileValue[1] & tmp$dominant_ctss.score<quantileValue[3] ),
    ]$dominant_ctss.pos),
  strand = tmp[which(tmp$dominant_ctss.score>quantileValue[1] & tmp$dominant_ctss.score<quantileValue[3] ),
  ]$strand,
  seqlengths = seqlengths(bsg))
  
  tmp3<-promoters(tmp3,upstream = 1,downstream = 1)
  tmp3<-tmp3[width(trim(tmp3)) == 2]
  
  
  c<-  melt(table(getSeq(bsg,tmp3)))
  #c<-c[c$Var1%notin%c("NN","NG"),]
  c<-c[!is.na(c$Var1),]
  d<-mutate(c,perc=value/sum(value))
  d$value<-NULL
  d<-d[order(d$perc,decreasing=F),]
  d$Var1<-factor(d$Var1,levels = d$Var1)
  
  ggplot(data=d,aes(x=factor(Var1),
                    y=(perc*100),fill=Var1))+
    geom_bar(stat="identity",colour = "black", size = 0.3,position=position_dodge())+
    ggtitle(paste0("summary of Dominant_TSS_2nd-3rdQuant_, ",sampleLabels(CAGEr_object)[i]," , -/+ 1 position"))+
    coord_flip()+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),text = element_text(size=20),legend.key.size = unit(0.8,"cm"))+
    labs(fill="")+
    scale_fill_viridis(option = "magma",discrete = TRUE)
  
  ggsave(filename = paste0("./Dominant_TSS_2to3Quant_",sampleLabels(CAGEr_object)[i],".pdf"),width = 12,height=10)
  
  
  tmp4<-GRanges(seqnames = tmp[which(tmp$dominant_ctss.score>quantileValue[3]),
  ]$seqnames,
  ranges =
    IRanges(start =  tmp[which(tmp$dominant_ctss.score>quantileValue[3] ),
    ]$dominant_ctss.pos,
    end =  tmp[which(tmp$dominant_ctss.score>quantileValue[3] ),
    ]$dominant_ctss.pos),
  strand =  tmp[which(tmp$dominant_ctss.score>quantileValue[3] ),
  ]$strand,
  seqlengths = seqlengths(bsg))
  
  tmp4<-promoters(tmp4,upstream = 1,downstream = 1)
  tmp4<-tmp4[width(trim(tmp4)) == 2]
  
  c<- melt(table(getSeq(bsg,tmp4)))
  #c<-c[c$Var1%notin%c("NN","NG"),]
  c<-c[!is.na(c$Var1),]
  
  d<-mutate(c,perc=value/sum(value))
  d$value<-NULL
  d<-d[order(d$perc,decreasing=F),]
  d$Var1<-factor(d$Var1,levels = d$Var1)
  
  ggplot(data=d,aes(x=factor(Var1),
                    y=(perc*100),fill=Var1))+
    geom_bar(stat="identity",colour = "black", size = 0.3,position=position_dodge())+
    ggtitle(paste0("summary of Dominant_TSS_4thQuant_, ",
                   sampleLabels(CAGEr_object)[i]," , -/+ 1 position"))+
    coord_flip()+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size=20),
          legend.key.size = unit(0.8,"cm"))+
    labs(fill="")+
    scale_fill_viridis(option = "magma",discrete = TRUE)

  ggsave(filename = paste0("./Dominant_TSS_4thQuant_,",sampleLabels(CAGEr_object)[i],".pdf"),width = 12,height=10)
  
  print(paste0(sampleLabels(CAGEr_object)[i],", plot: done"))
  }

 
 
 
 
 
 
 