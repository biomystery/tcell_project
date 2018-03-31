require(edgeR)
require(tidyverse)
require(ggplot2)
require(pheatmap)
pd <- readRDS(file = "./finalize/pd_fig.rds")
attach(pd)
# Fig.7A -maPlot ----------------------------------------------------------
par(mfrow=c(1,2))
for(i in 1:2){
  plotSmear(lrt[[i]],de.tags = rownames(y)[as.logical(edger.DEGs[[i]])],
            main=colnames(edger.comp.mat)[i],deCol="#90EE90")
  abline(h=c(-th.lfc, th.lfc), col="#90EE90")
  text(10,4,paste0('Up:',summary(edger.DEGs[[i]])[3],'\n Down:',summary(edger.DEGs[[i]])[1]))
}
par(mfrow=c(1,1))

# Fig.7B - venn diagam ----------------------------------------------------
require(venn)
venn(lapply(edger.DEGs,function(x) which(x==1)),zcolor = "style",ellipse = T,cexil = 0.75)


# Fig. 7C kegg_pathway analysis  ----------------------------------------------------------
set.seed(111)
bg_list <- substr(sample(rownames(genes.anno),length(glist.NES.anno$ensembl_gene_id)),1,18)
write.csv(file='./finalize/KEGG_pathway/bg_list.csv',bg_list,quote = F,row.names = F,col.names = NULL)
write.csv(file='./finalize/KEGG_pathway/nes_list.csv',glist.NES.anno$ensembl_gene_id,
          quote = F,row.names = F,col.names = NULL)
pd.kegg <- read.delim(file = './finalize/KEGG_pathway/david_go_nes_10_11_2017.txt')
pd.kegg <- pd.kegg%>%filter(Benjamini<0.05) %>% arrange(Fold.Enrichment)
pd.kegg$Term <- factor(pd.kegg$Term,levels = pd.kegg$Term)
p <- ggplot(pd.kegg,aes(y=Fold.Enrichment,x=Term))+ geom_bar(stat = 'identity',fill="#90EE90")+ coord_flip()+
  theme_bw()+theme(axis.title = element_blank(),text = element_text(size = 15))
ggsave(filename = "./finalize/subfigs/Subfig7c.eps",width = 3,height = 2.5,plot = p,
       scale = 2.5)

# Fig.7D heatmap  ---------------------------------------------------------
glist.NES.clust1 <- read.csv(file='./analysis3/08Homer/glist_clust1.txt',
                             stringsAsFactors = F,header = F)
glist.NES.clust2 <- read.csv(file='./analysis3/08Homer/glist_clust2.txt',
                             stringsAsFactors = F,header = F)
glist.NES.clust3 <- read.csv(file='./analysis3/08Homer/glist_clust3.txt',
                             stringsAsFactors = F,header = F)

pd.fig <- readRDS(file="./finalize/pd_fig.rds")
pd <- pd.fig$d_hm
setEPS()

h <- pheatmap(pd[,-ncol(pd)],scale = "row",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 6,#fontsize_col = 4,
              #color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              color= colorRampPalette(c("white", "#554266"))(10),#fontsize = 2,
              #treeheight_row = 1,treeheight_col = 10
              width = 2.6597,height = 3.125,
              breaks = seq(-m,m,length.out = 11),border_color = NA,#cellwidth = 2,
              gaps_col = 9,cutree_rows = 3)
ngene_cluster<- c(83,14,30)
rgap <- cumsum(ngene_cluster)
ord <- h$tree_row$order
pd.1 <- pd[ord,]
pd.1 <- pd.1[c((rgap[2]+1):rgap[3],(rgap[1]+1):rgap[2],1:rgap[1]),]
pd.2 <- t(scale(t(pd.1[,-ncol(pd.1)])))
scale.2 <- function(x) (x-min(x))/(max(x)-min(x))
if(F){for(i in 1:3){
  colns <- seq(i,18,by=3)
  pd.2[,colns] <- apply(pd.1[,colns],1,scale.2)
}
}
pd.2 <- t(apply(pd.1[,-ncol(pd.1)],1,scale.2))

if(T){setEPS()
postscript(file="./finalize/subfigs/subfig7d.eps",width = 2,height = 5)
range(pd.2)
pheatmap(pd.2,scale = "none",show_rownames = F,
              cluster_cols = F,cluster_rows = F,fontsize_row = 5,#fontsize_col = 4,
         breaks = seq(-1.73,2.82,length.out = 11),
              #color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              color= colorRampPalette(c("white", "#554266"))(10),#fontsize = 2,
              #treeheight_row = 1,treeheight_col = 10
              width = 2.6597*3,height = 3.125*3,show_colnames = F,legend = F,
              border_color = NA,cutree_rows = 3,#cellwidth = 2,
              gaps_col = 9,gaps_row = cumsum(ngene_cluster[c(3,2,1)]))

dev.off()
}
pd.fig$d_hm <- pd

# Fig.7EF-motifs ----------------------------------------------------------

# 1. logos
require(ggseqlogo)
require(ggplot2)
require(gridExtra)
attach(pd$ef_motifs_pwa)

p.c50 <- ggplot() + geom_logo( t(pd$ef_motifs_pwa[[1]]),method='prob' ) + theme_logo() + 
  theme(text = element_blank())

p.crel <- ggplot() + geom_logo( t(pd$ef_motifs_pwa[[2]]),method='prob' ) + theme_logo() + 
  theme(text = element_blank())
p.a50 <- ggplot() + geom_logo( t(pd$ef_motifs_pwa[[3]]),method='prob' ) + theme_logo() + 
  theme(text = element_blank())
if(T){
  w = 1;h=.25;s=2
  ggsave(filename = "./finalize/subfigs/Subfig7_a50.eps",width = w,height = h,plot = p.a50,scale = s)
  ggsave(filename = "./finalize/subfigs/Subfig7_c50.eps",width = w,height = h,plot = p.c50,scale = s)
  ggsave(filename = "./finalize/subfigs/Subfig7_crel.eps",width = w,height = h,plot = p.crel,scale = s)
}



# Fig.7EF-motifs-bar ----------------------------------------------------------
pd <- read.csv(file = "./analysis3/wholeMM10/all.anno.reduced.csv",stringsAsFactors = F)
#glist.NES <- read.csv(file='./analysis3/08Homer/glist_clustAll.txt',
#                      stringsAsFactors = F,header = F)

pd <- pd%>%mutate(gcate = ifelse(Nearest.Ensembl %in% glist.NES.clust1$V1,"Late",
                    ifelse(Nearest.Ensembl %in% glist.NES.clust2$V1,"Early & persistent",
                           ifelse(Nearest.Ensembl %in% glist.NES.clust3$V1,"Early","Background"))))%>%
  distinct()
if(T){
  th <-8
  pd <- pd %>%   mutate(p65_p50=p65_p50_mouse.Best.Motif.log.odds.Score>=th,
                        crel_p50=crel_p50_mouse.Best.Motif.log.odds.Score>=th,
                        crel=crel_mouse.Best.Motif.log.odds.Score>=th)
  
  
  pd.sum <- pd %>% group_by(gcate)%>%
    summarise(crel_p50 = sum(crel_p50)/n()*100,
              p65_p50 = sum(p65_p50)/n()*100,
              crel = sum(crel)/n()*100)%>% 
    gather(key = "TF",value = "percent",c(2:4))
  
  pd.sum$gcate <- factor(pd.sum$gcate,levels = c("Early", "Early & persistent", "Late" , "Background"))
  cols <- c(rep('black',3),"grey")
  ggplot(pd.sum, aes(x=gcate,y=percent,fill=gcate))+ geom_bar(stat = 'identity',width =.8)+
    scale_fill_manual(values = cols) + facet_wrap(~ TF) + ylim(0,40)

  sapply(unique(pd.sum$TF),function(x){
    p <- ggplot(pd.sum%>% filter(TF==x), aes(x=gcate,y=percent,fill=gcate))+ geom_bar(stat = 'identity',width =.8)+
      scale_fill_manual(values = cols) + ylim(0,40) + theme_bw() + theme(title = element_blank(),legend.position = "none",
                                                                         axis.text =  element_blank())
    ggsave(filename = paste0("./finalize/subfigs/subfig7_menrich",x,".eps"),p,
           width = 1,height = 1.2,scale = 2.5)
  })  
}

gcates <- c("Early", "Early & persistent", "Late" , "Background")
sapply(1:3,function(x){
  tb <- table(pd%>% filter(gcate %in% gcates[c(x,4)])%>%
                select(crel,gcate))
  fisher.test(tb,alternative="greater")
  })

 