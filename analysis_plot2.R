load('./.RData')
# Plot number of gene vs. threshold --------------------------------------------------------
# glist.wtup - fdr 0.5, lfc = 1 for wt LRT test 

load(file='./2017-04-04-analysis1/edgeR_cmp.data')

# compare pairs 
edger.comp.mat <- data.frame(
  NESvWT_12_0= c(1,-1,0,0,0,0),
  NESvWT_12_12= c(0,0,1,-1,0,0),
  NESvWT_24_24=c(0,0,0,0,1,-1))

lrt <- sapply(colnames(edger.comp.mat),function(x)
  glmLRT(fit,contrast = edger.comp.mat[,x]),simplify = F)


ngeneFunc<- function(th.fdr,th.lfc) {
  edger.NESvWT.DEGs <- lapply(lrt,function(x) 
    decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))

  # NESvWT.glist
  g.cates <- lapply(edger.NESvWT.DEGs,function(x)which(x==-1))
  NESvWT.glist <- data.frame( ids = unique(unlist(g.cates)))
  NESvWT.glist<- NESvWT.glist %>% mutate(ensemble = y$genes$Geneid[ids])
  rownames(NESvWT.glist) <- as.character(NESvWT.glist$ids)
  
  tmp <- intersect(NESvWT.glist$ids,glist.wtup)
  pd <- rpkm[NESvWT.glist[as.character(tmp),]$ensemble,]
  
  pd<- pd[apply(pd,1,max)>=2,]
}

pd.ngene.fdr05 <- data.frame(lfc= -seq(0,1,by=.1),
                             ngenes=sapply(seq(0,1,by=.1),function(x) nrow(ngeneFunc(0.05,th.lfc = x))))
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =4,hjust=1)



# plot WT & NES red genes (heatmap) ----------------------------------------------------------
# load tmp.dic from .RData 
#gene.dic <- data.frame(Symbol=getGsymble(y$genes$Geneid),
#                       stringsAsFactors = F)
#write.csv(file = './data/geneDic.csv',gene.dic,quote = F,row.names = F)
gene.dic <- read.csv(file="./data/geneDic.csv",stringsAsFactors = F,row.names = 1)
pd <- ngeneFunc(th.fdr = 0.05,th.lfc = 0.5)
pd$Ensembl <- rownames(pd)
rownames(pd) <- gene.dic[rownames(pd),1]
m <- max(abs(range(t(scale(t(pd[,-ncol(pd)]))))))
ord <- c(1:3,7:9,13:15,4:6,10:12,16:18)
pd <- pd[,c(ord,ncol(pd))]
if(T){
  # need load data first 
  # save data for the first time usage 
  setEPS()
  postscript(file="./finalize/subfig7d.eps",width = 2.6597,height = 3.125)
  h <- pheatmap(pd[,-ncol(pd)],scale = "row",show_rownames = T,
                cluster_cols = F,cluster_rows = T,fontsize_row = 6,#fontsize_col = 4,
                #color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
                color= colorRampPalette(c("white", "#554266"))(10),#fontsize = 2,
                #treeheight_row = 1,treeheight_col = 10
                width = 2.6597,height = 3.125,
                breaks = seq(-m,m,length.out = 11),border_color = NA,#cellwidth = 2,
                gaps_col = 9,cutree_rows = 3)
  dev.off()
  pd.fig <- readRDS(file="./finalize/pd_fig.rds")
  pd.fig$d_hm <- pd
  saveRDS(pd.fig,file = "./finalize/pd_fig.rds")
  
}


# export three list 
gs <- rownames(pd)[h$tree_row$order]
idx1 <- grep('Podnl1',gs)
idx2 <- grep('Ccl1',gs)
pd <- data.frame(pd,
                 clust=c(rep(1,idx1),rep(2,idx2-idx1),rep(3,length(gs)-idx2)))
write.csv(file = './data/rpkm_NESgene_fdr05_lfc.5.csv',pd)

sapply(1:3,function(i) write.table(pd$Ensembl[pd$clust==i],
                                   file = paste0('./data/glist_clust',i,'.txt'),
                                   sep = '\n',quote = F))

# plot WT & NES red genes (bar)  ----------------------------------------------------------
pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
rownames(pd) <- tmp.dic[as.character(glist.NESred.2),]$Symbol
gs <-names(which(apply(pd,1,max)>2))
gs <- gs[order(gs)]
barfunMerge(names(which(apply(pd,1,max)>2)))

pd.2 <- pd[which(apply(pd,1,max)>2),ord]
h <- pheatmap(pd[which(apply(pd,1,max)>2),ord],scale = "row",show_rownames = T,
         cluster_cols = F,cluster_rows = T,fontsize_row = 4,
         color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
         breaks = seq(-m,m,length.out = 21),border_color = NA,
         gaps_col = 9,cutree_rows = 2)
sep <- which(rownames(pd.2)[h$tree_row$order]=='Irf5')
barfunMerge.2 <- function(sep.1=1,sep.2=sep){
  p <- barfunMerge(rownames(pd.2)[h$tree_row$order[sep.1:sep.2]])  
  # adding the testing results 
  p.dat <- p$data
  tmp.func<-function(g=as.character(p.dat$gene[idx])){
    ids <- NESvWT.glist$ids[(NESvWT.glist$Symbol==g)]
    data.frame(treat = sub("NESvWT_","",names(NESvWT.glist.l)),
               cate = sapply(NESvWT.glist.l,function(x) ids %in% x),
               gene = g)
    
  }
  p.stat.dat.logic <- lapply(as.list(as.character(unique(p.dat$gene))),
                             tmp.func)
  p.stat.dat.logic <- do.call(rbind,p.stat.dat.logic)
  p.stat.dat.logic <- p.stat.dat.logic %>% arrange(treat)
  p.stat.dat <- p.dat%>% 
    group_by(treat,gene) %>%
    summarise(mm = max(m)+max(s)+1)%>% 
    arrange(treat,as.character(gene))
  
  # put NA if not significant 
  p.stat.dat$mm[!p.stat.dat.logic$cate] =NA
  p.stat.dat<- p.stat.dat%>% mutate(geno="NES")
  # added statistics into the plot 
  p + geom_point(data = p.stat.dat,
                 aes(treat,mm),shape='*',size=6,colour='red')
  
}
barfunMerge.2()
barfunMerge.2(sep.1 = sep+1,
              sep.2 = nrow(pd.2))

 