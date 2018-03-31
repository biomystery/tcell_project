getEntrezID <- function(glist){
  
  try(symbol <- mapIds(org.Mm.eg.db,
                       keys=glist,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first"))
  
  ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")
  if(!exists('symbol'))
    symbol <- getBM(attributes=c('ensembl_gene_id','entrezgene'),filters = 'ensembl_gene_id',
                    values  =glist, mart = ensembl)$entrezgene
  
  if( sum(is.na(symbol))>0)
    symbol[is.na(symbol)] <- getBM(attributes=c('ensembl_gene_id','entrezgene'),filters = 'ensembl_gene_id',
                                   values  =glist[is.na(symbol)], mart = ensembl)$entrezgene
  
  if(sum(is.na(symbol))>0)
    symbol[(is.na(symbol))]<- names(symbol[(is.na(symbol))])
  
  symbol
  
}

# seqlogo  ----------------------------------------------------------------
require(ggseqlogo)
require(gridExtra)

data(ggseqlogo_sample)
ggplot() + geom_logo( seqs_dna$sample_dna_1 ) + theme_logo()

colns <- c('A',"C","G","T")
colns.rev <- c("T","G","C","A")
readLines('./analysis3/08Homer/Siggers_PWMs/siggersknown.motif',12)
pwm.c50<- read.table(file='./analysis3/08Homer/Siggers_PWMs/siggersknown.motif',skip = 11,
           nrows = 11,col.names = colns)
pwm.cRel <- read.table(file='./analysis3/08Homer/Siggers_PWMs/siggersknown.motif',skip = 1,
                                 nrows = 9,col.names = colns)
pwm.a50 <- read.table(file='./analysis3/08Homer/Siggers_PWMs/siggersknown.motif',skip = 33,
                       nrows = 12,col.names = colns)
pwm.c50.rev <- pwm.c50[rev(rownames(pwm.c50)),]
colnames(pwm.c50.rev)<- colns.rev
p.c50 <- ggplot() + geom_logo( t(pwm.c50),method='prob' ) + theme_logo()
p.a50 <- ggplot() + geom_logo( t(pwm.a50),method='prob' ) + theme_logo()
p.crel <- ggplot() + geom_logo( t(pwm.cRel),method='prob' ) + theme_logo()

pd.fig$ef_motifs_pwa <- list(pwm.c50,pwm.cRel,pwm.a50)

p.c50.rev <- ggplot() + geom_logo( t(pwm.c50.rev),method='prob' ) + theme_logo()
grid.arrange(p.a50,p.c50,p.crel,ncol=1)

# creat three gene list  --------------------------------------------------

# write glist for three clusts 
pd <- read.csv(file = './data/rpkm_NESgene_fdr05_lfc.5.csv',
               stringsAsFactors = F,row.names = 1)

sapply(1:3,function(i) write.table(pd$Ensembl[pd$clust==i],
                                   file = paste0('./data/glist_clust',i,'.txt'),
                                   sep = '\n',quote = F,row.names = F,col.names =F))

# run homer @ SeqC 

# homer call2: all clust genes vs remaining induced genes for knownmotif 
# write gene list for background (WT induced genes - NES dependent glist)
write.table(file='./data/WT_induced_remain_genes.txt',
            setdiff(y$genes$Geneid[glist.wtup],pd$Ensembl),
            sep = '\n',quote = F,row.names = F,col.names = F)

# filtering out the rpkm <=2 
tmp.rpkm <- rpkm[y$gene$Geneid[glist.wtup],]
tmp.rpkm<- subset(tmp.rpkm,apply(tmp.rpkm,1,function(x) max(x)>2)) # 1012 genes of 1597 left 
tmp.rpkm$Ensembl <- rownames(tmp.rpkm)
write.csv(file='./data/rpkm_wt_induced_FDR05_LFC1_filtering_rpkm2.txt',tmp.rpkm)
all(pd$Ensembl %in%rownames(tmp.rpkm))
write.table(file='./data/WT_induced_remain_genes_filtering_rpkm2.txt',
            setdiff(rownames(tmp.rpkm),pd$Ensembl),
            sep = '\n',quote = F,row.names = F,col.names = F)


# homer call 3: for each individual clust call against the remaining genes in all list 
sapply(1:3,function(i) write.table(file = paste0('./data/glist_clust',i,'_bg.txt'),
                                   setdiff(rownames(tmp.rpkm),pd$Ensembl[pd$clust==i]),
                                   sep = '\n',quote = F,row.names = F,col.names = F))




#(findmotifs.pl) plot motif enrichment result-------------------------------------------
require(ggplot2);require(dplyr); require(tidyr)
motif.res <- read.table(file = './analysis3/08Homer/glist_clustAll/knownResults.txt',
                        sep = '\t',skip = 1)

motif.res <- motif.res[,c(1,4,7,9)]
motif.res <- motif.res%>% gather(key = 'type',value = 'percent',3:4)
motif.res<- motif.res%>% mutate(percent=as.numeric(sub('%','',motif.res$percent)))
ggplot(motif.res%>% filter(V1 %in% c('crel_mouse','p65_mouse')),aes(V1,percent)) + 
  geom_bar(stat = 'identity',aes(fill=type),position = 'dodge',colour='white')

# annotatePeaks motif res -------------------------------------------------

motifs.c50 <- read.csv(file = './analysis3/wholeMM10/c50.anno.reduced.csv',
                       stringsAsFactors = F)
motifs.crel <- read.csv(file = "./analysis3/wholeMM10/")
motifs.a50 <- read.csv(file = './analysis3/wholeMM10/a50.anno.reduced.csv',
                       stringsAsFactors = F)
all(motifs.a50$Gene.Name %in% motifs.c50$Gene.Name) # TRUE 

motifs.a50<- motifs.a50 %>%
  arrange(Gene.Name,desc(p65_p50_mouse.Best.Motif.log.odds.Score)) %>%
  distinct(Gene.Name,.keep_all = T) # reomve duplicate and keep best

motifs.c50<- motifs.c50 %>%
  arrange(Gene.Name,desc(crel_p50_mouse.Best.Motif.log.odds.Score)) %>%
  distinct(Gene.Name,.keep_all = T) # reomve duplicate and keep best

motifs.genome <- (right_join(motifs.a50,motifs.c50[,c(2,4)],"Gene.Name"))

 # whole genome plotting  ------------------------------------------------


motifs.genome.2 <- motifs.genome #backup 
motifs.genome$p65_p50_mouse.Best.Motif.log.odds.Score<-motifs.genome$p65_p50_mouse.Best.Motif.log.odds.Score/12
motifs.genome$crel_p50_mouse.Best.Motif.log.odds.Score<-motifs.genome$crel_p50_mouse.Best.Motif.log.odds.Score/11

if(T){
pm<-ggplot(motifs.genome,aes(p65_p50_mouse.Best.Motif.log.odds.Score,
                         crel_p50_mouse.Best.Motif.log.odds.Score))+
  stat_density2d(aes(fill = ..density..), 
                 geom = "tile", contour = FALSE, n = 100) +
  geom_point(alpha=0.2,size=1)+
  geom_point(data = subset(motifs.genome,
                           isNES==4),aes(colour=isNES)) +
  geom_abline(slope = 1,intercept = 0)+
  scale_fill_continuous(low = "white", high = "red")+
  geom_smooth(method = "lm")+
  theme_bw()+theme(legend.position = 'none')

ptop <- ggplot(motifs.genome,aes(p65_p50_mouse.Best.Motif.log.odds.Score)) + 
  geom_density(alpha=0.2)+ theme_bw()+
  theme(axis.title.x =element_blank())

pleft<-ggplot(motifs.genome,aes(crel_p50_mouse.Best.Motif.log.odds.Score)) + 
  geom_density(alpha=0.2)+theme_bw()+
  coord_flip()+theme(legend.position = 'none')+
  scale_y_reverse()

pempty <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
require(gridExtra)
p <- grid.arrange(pempty,ptop , pleft,pm, 
                  ncol=2, nrow=2, widths=c(1.4, 4), heights=c(1.4, 4))
}

motifs.genome <- motifs.genome.2


# motif: WT list vs. glist NES  -----------------------------------------
glist.WT <- read.csv(file='./analysis3/08Homer/WT_induced_genes_rpkm2.txt',
                     stringsAsFactors = F,header = F)

motifs.genome <- subset(motifs.genome,Nearest.Ensembl%in% glist.WT$V1)
dim(motifs.genome) #930 

glist.NES <- read.csv(file='./analysis3/08Homer/glist_clustAll.txt',
                      stringsAsFactors = F,header = F)
motifs.genome<-motifs.genome%>% mutate(isNES =Nearest.Ensembl%in%glist.NES$V1)#110

if(T){
  pm<-ggplot(motifs.genome,aes(p65_p50_mouse.Best.Motif.log.odds.Score,
                               crel_p50_mouse.Best.Motif.log.odds.Score))+
    #stat_density2d(aes(fill = ..density..), 
    #               geom = "tile", contour = FALSE, n = 100) +
    geom_point(alpha=0.4,aes(colour=isNES))+
    geom_abline(slope = 1,intercept = 0)+
    scale_fill_continuous(low = "white", high = "red")+
    geom_smooth(method = "lm")+
    theme_bw()+theme()
  
  ptop <- ggplot(motifs.genome,aes(p65_p50_mouse.Best.Motif.log.odds.Score)) + 
    geom_density(alpha=0.2,aes(fill=isNES))+ theme_bw()+
    theme(axis.title.x =element_blank())
  
  pleft<-ggplot(motifs.genome,aes(crel_p50_mouse.Best.Motif.log.odds.Score)) + 
    geom_density(alpha=0.2,aes(fill=isNES))+theme_bw()+
    coord_flip()+theme(legend.position = 'none')+
    scale_y_reverse()
  
  pempty <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  require(gridExtra)
  p <- grid.arrange(pempty,ptop , pleft,pm, 
                    ncol=2, nrow=2, widths=c(1.4, 4), heights=c(1.4, 4))
}

# boxplot for a50 and c50 score changes   ---------------------------------

motifs.genome.sum <- (motifs.genome%>% gather(key=TF,value=mlogOddsSore,4:5))
ggplot(motifs.genome.sum,aes(isNES,mlogOddsSore)) + geom_boxplot(aes(fill=isNES))+
  facet_wrap(~TF)

x=with(motifs.genome.sum,
     mlogOddsSore[TF=="p65_p50_mouse.Best.Motif.log.odds.Score"&isNES])
y=with(motifs.genome.sum,
     mlogOddsSore[TF=="p65_p50_mouse.Best.Motif.log.odds.Score"&!isNES])
a50.test <- t.test(x,y,
       alternative='greater') #p-value = 0.0001317
a50.kstest <- ks.test(x,y,
                   alternative='greater') #p-value = 0.0001317


x=with(motifs.genome.sum,
       mlogOddsSore[TF=="crel_p50_mouse.Best.Motif.log.odds.Score"&isNES])
y=with(motifs.genome.sum,
       mlogOddsSore[TF=="crel_p50_mouse.Best.Motif.log.odds.Score"&!isNES])
c50.test <- t.test(x,y,
                   alternative='greater') #p-value = 0.0001317
c50.kstest <- ks.test(x,y,
                   alternative='greater') #p-value = 0.0001317


# compare with background  ------------------------------------------------

motifs.genome.2 <- motifs.genome.2 %>% mutate(gcate=ifelse(Nearest.Ensembl %in% glist.NES$V1,'NES gene',
                                        ifelse(Nearest.Ensembl %in% glist.WT$V1,'WT gene-NES gene',
                                                      'Background')))
motifs.genome.2$gcate <- factor(motifs.genome.2$gcate,c('WT gene-NES gene',"NES gene",'Background'))
require(scales)
cols <- c(hue_pal()(2),'black')
p1 <- ggplot(motifs.genome.2,aes(crel_p50_mouse.Best.Motif.log.odds.Score)) +
  geom_density(aes(fill=gcate),alpha=0.4) + 
  scale_fill_manual(values = cols)
p2 <- ggplot(motifs.genome.2,aes(p65_p50_mouse.Best.Motif.log.odds.Score)) + 
  geom_density(aes(fill=gcate),alpha=0.4)+ 
  scale_fill_manual(values = cols)
grid.arrange(p1,p2,ncol=1)

#update to individual clusters 
glist.NES.clust1 <- read.csv(file='./analysis3/08Homer/glist_clust1.txt',
                             stringsAsFactors = F,header = F)
glist.NES.clust2 <- read.csv(file='./analysis3/08Homer/glist_clust2.txt',
                             stringsAsFactors = F,header = F)
glist.NES.clust3 <- read.csv(file='./analysis3/08Homer/glist_clust3.txt',
                             stringsAsFactors = F,header = F)


motifs.genome.3 <- motifs.genome.2
motifs.genome.3$gcate <- as.character(motifs.genome.3$gcate)
motifs.genome.3$gcate[motifs.genome.3$Nearest.Ensembl %in% glist.NES.clust1$V1]='NES late'  
motifs.genome.3$gcate[motifs.genome.3$Nearest.Ensembl %in% glist.NES.clust2$V1]='NES early and transient'  
motifs.genome.3$gcate[motifs.genome.3$Nearest.Ensembl %in% glist.NES.clust3$V1]='NES early'  

cols <- c('black',hue_pal()(4))
p1 <- ggplot(motifs.genome.3 
        ,aes(crel_p50_mouse.Best.Motif.log.odds.Score)) +
  geom_density(aes(fill=gcate),alpha=0.2)+ 
  scale_fill_manual(values = cols)
p2 <- ggplot(motifs.genome.3 
             ,aes(p65_p50_mouse.Best.Motif.log.odds.Score)) +
  geom_density(aes(fill=gcate),alpha=0.2)+ 
  scale_fill_manual(values = cols)
  
grid.arrange(p1,p2,ncol=2)
# enrichment calculation --------------------------------------------------
# examples 
bg <- 871; bg.t <- 7; target <- 110;target.t <- 6
bg <- 22322; bg.t <- 55; target <- 110;target.t <- 8

fisher.test(matrix(c(target.t,target,bg.t,bg),
                   nrow=2,ncol=2),alternative="greater")
motifs.genome.2<- motifs.genome.2%>% 
  mutate(p65_p50=p65_p50_mouse.Best.Motif.log.odds.Score>=8,
         crel_p50=crel_p50_mouse.Best.Motif.log.odds.Score>=8)

tb <- table(motifs.genome.2%>% filter(gcate!='NES gene')%>%
  select(p65_p50,gcate)%>%
    mutate(gcate=droplevels(gcate)))
fisher.test(tb,alternative="greater")

tb <- table(motifs.genome.2%>% filter(gcate!='WT gene-NES gene')%>%
              select(p65_p50,gcate)%>%
              mutate(gcate=droplevels(gcate)))
fisher.test(tb[rev(rownames(tb)),],alternative="greater")

tb <- table(motifs.genome.2%>% filter(gcate!='NES gene')%>%
              select(crel_p50,gcate)%>%
              mutate(gcate=droplevels(gcate)))
fisher.test(tb,alternative="greater")

tb <- table(motifs.genome.2%>% filter(gcate!='WT gene-NES gene')%>%
              select(crel_p50,gcate)%>%
              mutate(gcate=droplevels(gcate)))
fisher.test(tb[rev(rownames(tb)),],alternative="greater")

tb <- table(motifs.genome.2%>% filter(gcate=='NES gene')%>%
              select(p65_p50,crel_p50,gcate)%>%
              mutate(gcate=droplevels(gcate)))
fisher.test(tb[rev(rownames(tb))],alternative="greater")

# overlay 
pd <- sapply(c('NES late',"NES early and transient",
               "NES early"),function(x)
                 nrow(motifs.genome.2%>% filter(gcate==x & crel_p50)))
pie(pd,labels = paste(names(pd),round(pd/sum(pd)*100),"%"))

# for each class
cols <- c('grey',rep('black',4))
ggplot(motifs.genome.3%>% group_by(gcate)%>%
  summarise(crel_p50 = sum(crel_p50)/n()*100,
            p65_p50 = sum(p65_p50)/n()*100)%>% 
  gather(key = "TF",value = "percent",c(2:3)),
  aes(gcate,percent,fill=gcate))+ geom_bar(stat = 'identity')+
  facet_wrap(~TF)+ scale_fill_manual(values = cols)


# bar plot 
tb <- table(motifs.genome.2%>% 
              select(crel_p50,gcate)%>%
              mutate(gcate=droplevels(gcate)))
tb <- as.data.frame(tb)

tb.2 <- table(motifs.genome.2%>% 
              select(p65_p50,gcate)%>%
              mutate(gcate=droplevels(gcate)))
tb.2 <- as.data.frame(tb.2)

pd.fig <- readRDS(file = "./finalize/pd_fig.rds")

tb<- rbind(data.frame(tb[,-1],TF='crel:p50'),
      data.frame(tb.2[,-1],TF='p65:p50'))

ggplot(tb %>% group_by(gcate,TF) %>%
  summarise(percent=last(Freq)/sum(Freq)*100),
  aes(gcate,percent,fill=gcate))+ geom_bar(stat = 'identity')+
  scale_fill_manual(values = cols)+facet_wrap(~TF)

pd.fig$ef_motifs <- tb
saveRDS(file="./finalize/pd_fig.rds",object = pd.fig)
# plot GOterm for glist  --------------------------------------------------
gopro <- read.table(file="./analysis3/08Homer/glist_clustAll/biological_process2.txt",
                    sep = '\t',nrows = 15,quote = T)
gopro$V1 <- factor(gopro$v1,levels = gopro$V1)
barplot(gopro$V4)
ggplot(gopro,aes(V1,V4))+ geom_bar(stat = 'identity')
write.csv(gopro,file='headgo.csv',quote = F)

gopro<-read.csv(file='headgo.csv',quote = F,stringsAsFactors = F)

# panther -----------------------------------------------------------------
readLines('./analysis3/Panther/analysis.txt',12)
go.res<- read.table(file='./analysis3/Panther/analysis.txt',skip = 11,
            header = T,sep = '\t')

require(ggplot2)
require(RColorBrewer)
go.res$glist_clustAll.txt..fold.Enrichment. <- as.numeric(as.character(go.res$glist_clustAll.txt..fold.Enrichment.))
go.res$glist_clustAll.txt..fold.Enrichment.[1] <-.2
bks <- c(-0.0001,1.30103,2,3,4,5,8)
go.res$pcut <-   cut(-log10(go.res$glist_clustAll.txt..P.value.),bks)
cols <- colorRampPalette(c('white','red4'))(length(bks)-1)[2:(length(bks)-1)]

ggplot(go.res,aes(reorder(PANTHER.GO.Slim.Biological.Process,glist_clustAll.txt..fold.Enrichment.),glist_clustAll.txt..fold.Enrichment.))+
  geom_bar(stat = 'identity',aes(fill=-log10(glist_clustAll.txt..P.value.)))+theme_bw() +coord_flip()+
  theme(axis.title.y = element_blank())+
  scale_fill_continuous('-log10(pvalue)',high = 'red',low = 'grey80')+
  ylab("Fold enrichment") + 
  geom_hline(yintercept = 1,linetype=2,colour=grey(.8))


# WT and NES list against whole genome 
go.res<- read.table(file='./analysis3/Panther/WT_NES_list_whole_ref.txt',skip = 11,
                    header = T,sep = '\t')
colnames(go.res)
attach(go.res)
go.res$WT.list..fold.Enrichment.<- as.numeric(as.character(go.res$WT.list..fold.Enrichment.))
go.res$NES.list..fold.Enrichment.<- as.numeric(as.character(go.res$NES.list..fold.Enrichment.))
go.res[is.na(go.res)]<- 0.2
go.res.sum<- go.res%>% select(c(1,6,11))%>% 
  gather(key = "glist",value="FoldEnrichment",2:3)

go.res.sum<- go.res%>% select(c(1,6,7,11,12))%>% 
 gather(key = "glist",value="FoldEnrichment",2:4)

lvs <- (go.res%>% arrange(WT.list..fold.Enrichment.))$PANTHER.GO.Slim.Biological.Process
  
go.res.sum$PANTHER.GO.Slim.Biological.Process <- factor(go.res.sum$PANTHER.GO.Slim.Biological.Process,
                                                        levels = lvs)

ggplot(go.res.sum,aes(PANTHER.GO.Slim.Biological.Process,FoldEnrichment))+
  geom_bar(aes(fill=glist),position = "dodge",stat = 'identity')+
  coord_flip()+ geom_hline(yintercept = 1,linetype=2,colour=grey(.8))+
  ylab("FoldEnrichment against whole genome")


ggplot(go.res.sum%>% filter(glist=="NES.list..fold.Enrichment.")%>%arrange(desc(FoldEnrichment)),
       aes(PANTHER.GO.Slim.Biological.Process,FoldEnrichment))+
  geom_bar(aes(fill=glist),position = "dodge",stat = 'identity')+
  coord_flip()+ geom_hline(yintercept = 1,linetype=2,colour=grey(.8))+
  ylab("FoldEnrichment against whole genome")

# heatmap of pvalue 
bks <- c(-0.0001,1.30103,2,3,4,5,8,10,30)
go.res.sum.2<- go.res%>% select(c(1,7,12))%>% 
  gather(key = "glist",value="pvalue",c(2,3))%>%
  mutate(pcut= cut(-log10(pvalue),breaks = bks))
lvs <- (go.res%>%arrange(desc(WT.list..P.value.)))$PANTHER.GO.Slim.Biological.Process
go.res.sum.2$PANTHER.GO.Slim.Biological.Process <- factor(go.res.sum.2$PANTHER.GO.Slim.Biological.Process,
                                                          levels = lvs)
ggplot(go.res.sum.2,aes(glist,PANTHER.GO.Slim.Biological.Process))+
  geom_tile(aes(fill=pcut))+
  scale_fill_manual("-log10(pvalue)",
                    values = colorRampPalette(c('white','red4'))(length(bks)-1))+
  theme_bw() +theme(axis.text.x = element_blank(),
                    axis.title=element_blank())

# export gene list  -------------------------------------------------------
glist.WT<-glist.WT%>%  mutate(gcate=ifelse(V1 %in% glist.NES$V1,'NES gene',
                                 ifelse(V1 %in% glist.WT$V1,'WT gene-NES gene',
                                        'Background')))

glist.WT$gcate[glist.WT$V1 %in% glist.NES.clust1$V1]='NES late'  
glist.WT$gcate[glist.WT$V1 %in% glist.NES.clust2$V1]='NES early and transient'  
glist.WT$gcate[glist.WT$V1 %in% glist.NES.clust3$V1]='NES early'  
summary(as.factor(glist.WT$gcate))


require(biomaRt)
mart <- useMart(biomart = 'ensembl',dataset = 'mmusculus_gene_ensembl')

gene_annotation_func<- function(glist){
  anno.list <- c("ensembl_gene_id","mgi_symbol","gene_biotype","description")
  getBM(attributes = anno.list,
        filters = anno.list[1],
        values = glist,
        mart = mart)
}

glist.WT.anno <- gene_annotation_func(glist.WT$V1)
rownames(glist.WT.anno) <- glist.WT.anno$ensembl_gene_id
rownames(glist.WT) <- glist.WT$V1
glist.WT.anno$gcate <- glist.WT[glist.WT.anno$ensembl_gene_id,'gcate']

for(i in 1:3) glist.WT.anno<-rbind(glist.WT.anno,NA)
glist.WT.anno$ensembl_gene_id[(nrow(glist.WT.anno)-2):nrow(glist.WT.anno)]<-
  setdiff(glist.WT$V1,glist.WT.anno$ensembl_gene_id)

write.csv(file="final.genelist.anno.csv",
          glist.WT.anno,row.names = F)


# final code 
write.csv(file="motif.res.csv",
          (motifs.genome.3%>% filter(gcate!='Background')))

  
