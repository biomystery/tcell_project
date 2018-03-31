head(p$data)
dim(p$data)
p.dat <- p$data
head(p.dat)
g <- p.dat$gene[1]
g
g <- as.character(p.dat$gene[1])
head(NESvWT.glist)
which(NESvWT.glist$Symbol==g)
NESvWT.glist$ids[(NESvWT.glist$Symbol==g)]
ids <- NESvWT.glist$ids[(NESvWT.glist$Symbol==g)]
lapply(NESvWT.glist.l,function(x) ids %in% x)
tmp.func <- function(x) ids %in% x

do.call(tmp.func,NESvWT.glist)
ids
lapply(NESvWT.glist.l,tmp.func)
unlist(lapply(NESvWT.glist.l,tmp.func))
sapply(NESvWT.glist.l,tmp.func)
data.frame(sapply(NESvWT.glist.l,tmp.func))
data.frame(g=sapply(NESvWT.glist.l,tmp.func))
data.frame(t(sapply(NESvWT.glist.l,tmp.func)))
p.dat%>% filter(gene==g)
p.dat%>% filter(gene==g) %>%
  group_by(treat) %>%
  summarise(mm = max(m))
p.dat%>%
  group_by(treat,g) %>%
  summarise(mm = max(m))
p.dat%>%
  group_by(treat,gene) %>%
  summarise(mm = max(m))
p.stat.dat <- p.dat%>%
  group_by(treat,gene) %>%
  summarise(mm = max(m))
head(p.stat.dat)
data.frame(t(sapply(NESvWT.glist.l,tmp.func)))
data.frame(sapply(NESvWT.glist.l,tmp.func))
data.frame(treat = sub("NESvWT_","",names(NESvWT.glist.l)),
           cate = sapply(NESvWT.glist.l,tmp.func),
           gene = g)
tmp.func<-function(idx=1){
  g <- as.character(p.dat$gene[idx])
  ids <- NESvWT.glist$ids[(NESvWT.glist$Symbol==g)]
  data.frame(treat = sub("NESvWT_","",names(NESvWT.glist.l)),
             cate = sapply(NESvWT.glist.l,function(x) ids %in% x),
             gene = g)
}
dim(p.dat)
tmp.func<-function(g=as.character(p.dat$gene[idx])){
  ids <- NESvWT.glist$ids[(NESvWT.glist$Symbol==g)]
  data.frame(treat = sub("NESvWT_","",names(NESvWT.glist.l)),
             cate = sapply(NESvWT.glist.l,function(x) ids %in% x),
             gene = g)
}
unique(p.dat$gene)
do.call(tmp.func,as.character(unique(p.dat$gene)))
as.list(unique(p.dat$gene))
do.call(tmp.func,as.list(as.character(unique(p.dat$gene))))
?do.call
do.call("tmp.func",as.list(as.character(unique(p.dat$gene))))
as.list(as.character(unique(p.dat$gene)))
tmp.func((as.character(unique(p.dat$gene)))[1])
tmp.func((as.character(unique(p.dat$gene)))[2])
p.stat.dat.logic <- lapply(as.list(as.character(unique(p.dat$gene))),
                           tmp.func)
p.stat.dat.logic <- do.call(rbind,p.stat.dat.logic)
head(p.stat.dat.logic)
head(p.stat.dat)
head(p.stat.dat.logic%>%arrange(treat))
head(p.stat.dat%>%arrange(treat,gene))
head(p.stat.dat%>%arrange(treat,as.character(gene)))
p.stat.dat <- p.dat%>%
  group_by(treat,gene) %>%
  summarise(mm = max(m))%>% arrange(treat,as.character(gene))
p.stat.dat.logic <- p.stat.dat.logic %>% arrange(treat)
head(p.stat.dat.logic)
tail(p.stat.dat.logic)
tail(p.stat.dat)
p.stat.dat$mm[!p.stat.dat.logic$cate] =NA
p.stat.dat <- p.dat%>%
  group_by(treat,gene) %>%
  summarise(mm = max(m)+1)%>%
  arrange(treat,as.character(gene))
p.stat.dat$mm[!p.stat.dat.logic$cate] =NA
p + geom_point(data = p.stat.dat,
               aes(treat,max(m)),shape='*',size=6,colour='cyan')
head(p.dat)
p.stat.dat%>% mutate(geno="NES")
p.stat.dat<- p.stat.dat%>% mutate(geno="NES")
p + geom_point(data = p.stat.dat,
               aes(treat,max(m)),shape='*',size=6,colour='cyan')
head(p.stat.dat)
p + geom_point(data = p.stat.dat,
               aes(treat,mm),shape='*',size=6,colour='cyan')
p + geom_point(data = p.stat.dat,
               aes(treat,mm),shape='*',size=6,colour='red')
head(p.dat)
p.stat.dat <- p.dat%>%
  group_by(treat,gene) %>%
  summarise(mm = max(m)+max(s)+1)%>%
  arrange(treat,as.character(gene))
p.stat.dat$mm[!p.stat.dat.logic$cate] =NA
p.stat.dat<- p.stat.dat%>% mutate(geno="NES")
p + geom_point(data = p.stat.dat,
               aes(treat,mm),shape='*',size=6,colour='red')
p <- barfunMerge(rownames(pd.2)[h$tree_row$order[1:27]])
p.dat <- p$data
p.stat.dat.logic <- lapply(as.list(as.character(unique(p.dat$gene))),
                           tmp.func)
p.stat.dat.logic <- do.call(rbind,p.stat.dat.logic)
p.stat.dat.logic <- p.stat.dat.logic %>% arrange(treat)
p.stat.dat <- p.dat%>%
  group_by(treat,gene) %>%
  summarise(mm = max(m)+max(s)+1)%>%
  arrange(treat,as.character(gene))
p.stat.dat$mm[!p.stat.dat.logic$cate] =NA
p.stat.dat<- p.stat.dat%>% mutate(geno="NES")
p + geom_point(data = p.stat.dat,
               aes(treat,mm),shape='*',size=6,colour='red')
th.fdr <- .05; th.lfc <- .5
edger.NESvWT.DEGs <- lapply(lrt,function(x) decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
glist.NESred <- unique(unlist(lapply(edger.NESvWT.DEGs,function(x) which(x==-1))))
par(mfrow=c(2,2))
for(i in 1:3){
  plotSmear(lrt[[i]],de.tags = rownames(y)[as.logical(edger.NESvWT.DEGs[[i]])],
            main=colnames(edger.comp.mat)[i])
  abline(h=c(-th.lfc, th.lfc), col="blue")
  text(10,3,paste0('Up:',summary(edger.NESvWT.DEGs[[i]])[3],'\n Down:',summary(edger.NESvWT.DEGs[[i]])[1]))
}
par(mfrow=c(1,1))
edger.comp.mat <- data.frame(
  NESvWT_12_0= c(1,-1,0,0,0,0),
  NESvWT_12_12= c(0,0,1,-1,0,0),
  NESvWT_24_24=c(0,0,0,0,1,-1))
lrt <- sapply(colnames(edger.comp.mat),function(x)
  glmLRT(fit,contrast = edger.comp.mat[,x]))
glt  <- sapply(colnames(edger.comp.mat),function(x)
  glmTreat(fit,contrast = edger.comp.mat[,x],lfc = 1))
th.fdr <- .05; th.lfc <- .5
edger.NESvWT.DEGs <- lapply(lrt,function(x) decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
glist.NESred <- unique(unlist(lapply(edger.NESvWT.DEGs,function(x) which(x==-1))))
par(mfrow=c(2,2))
for(i in 1:3){
  plotSmear(lrt[[i]],de.tags = rownames(y)[as.logical(edger.NESvWT.DEGs[[i]])],
            main=colnames(edger.comp.mat)[i])
  abline(h=c(-th.lfc, th.lfc), col="blue")
  text(10,3,paste0('Up:',summary(edger.NESvWT.DEGs[[i]])[3],'\n Down:',summary(edger.NESvWT.DEGs[[i]])[1]))
}
par(mfrow=c(1,1))
require(venn)
venn(lapply(edger.NESvWT.DEGs,function(x) which(x==1)),zcolor = "style",cexil = 0.75)
length(unique(unlist(lapply(edger.NESvWT.DEGs,function(x)which(x==1)))))
NESvWT.glist.l <- lapply(edger.NESvWT.DEGs,function(x)which(x==-1))
g.cates <- calculate.overlap(lapply(edger.NESvWT.DEGs,function(x)which(x==-1)))
require(VennDiagram)
g.cates <- calculate.overlap(lapply(edger.NESvWT.DEGs,function(x)which(x==-1)))
NESvWT.glist <- data.frame(
  ids = unlist(g.cates),
  cates = unlist(sapply(1:7,function(x) rep(x,length(g.cates[[x]]))))
)
NESvWT.glist
NESvWT.glist<- NESvWT.glist %>% mutate(ensemble = y$genes$Geneid[ids])
NESvWT.glist<-NESvWT.glist%>% mutate(Symbol=getGsymble(ensemble))
require(biomaRt)
require(biomaRt)
getGsymble <- function(glist){
  try(symbol <- mapIds(org.Mm.eg.db,
                       keys=glist,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first"))
  ensembl = useEnsembl(biomart="ensembl",dataset = "mmusculus_gene_ensembl")
  if(!exists('symbol'))
    symbol <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',
                    values  =glist, mart = ensembl)$mgi_symbol
  if( sum(is.na(symbol))>0)
    symbol[is.na(symbol)] <- getBM(attributes=c('ensembl_gene_id','mgi_symbol'),filters = 'ensembl_gene_id',
                                   values  =glist[is.na(symbol)], mart = ensembl)$mgi_symbol
  if(sum(is.na(symbol))>0)
    symbol[(is.na(symbol))]<- names(symbol[(is.na(symbol))])
  symbol
}
NESvWT.glist<-NESvWT.glist%>% mutate(Symbol=getGsymble(ensemble))
NESvWT.glist
require(org.Mm.eg.db)
NESvWT.glist<-NESvWT.glist%>% mutate(Symbol=getGsymble(ensemble))
NESvWT.glist
pd <- log2(rpkm[NESvWT.glist$ensemble,]+1)
require(pheatmap)
pd <- t(scale(t(pd))); m <- max(abs(range(pd)))
pheatmap(pd,scale = "none",show_rownames = F,cluster_cols = F,
         colorRampPalette(c("navy", "white", "firebrick3"))(20),
         breaks = seq(-m,m,length.out = 21))

pheatmap(pd,scale = "none",show_rownames = F,cluster_cols = F,
         colorRampPalette(c("navy", "white", "firebrick3"))(20),
         breaks = seq(-m,m,length.out = 21))
dim(pd)
dim(NESvWT.glist)
venn(lapply(edger.NESvWT.DEGs,function(x) which(x==-1)),zcolor = "style",cexil = 0.75)
venn(lapply(edger.NESvWT.DEGs,function(x) which(x==-1)),zcolor = "style",cexil = 0.75)

length(glist.wtup)
data.frame (NESred=length(glist.NESred),
            WTup = length(glist.wtup),
            NESred_WTup = sum(glist.NESred%in% glist.wtup))
NESvWT.glist.l2 <- lapply(NESvWT.glist.l,function(x) x[x%in%glist.wtup])
venn(NESvWT.glist.l2,zcolor = "style",cexil = 0.75)
68+33+15+30+21
venn(NESvWT.glist.l2,zcolor = "style",cexil = 0.75)
names(edger.DEGs)
glist.NESred.2 <- glist.NESred[glist.NESred%in% glist.wtup]
length(glist.NESred.2)
tmp.dic <- NESvWT.glist; rownames(tmp.dic) <- as.character(tmp.dic$ids)
pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
rownames(pd) <- tmp.dic[as.character(glist.NESred.2),]$Symbol
pd <- log2(pd+1)
pd <- t(scale(t(pd))); m <- max(abs(range(pd)))
ord <- c(1:3,7:9,13:15,4:6,10:12,16:18)
h <- pheatmap(pd[,ord],scale = "none",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 6,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9)
pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
rownames(pd) <- tmp.dic[as.character(glist.NESred.2),]$Symbol
barfunMerge(rownames(pd)[h$tree_row$order[14:16]])
gs <-names(which(apply(pd,1,max)>2))
gs <- gs[order(gs)]
length(gs)
head(pd)
pd.2 <- pd[which(apply(pd,1,max)>2),ord]
h <- pheatmap(pd[which(apply(pd,1,max)>2),ord],scale = "row",show_rownames = T,
              cluster_cols = F,cluster_rows = T,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9,cutree_rows = 2)
h <- pheatmap(pd[which(apply(pd,1,max)>2),ord],scale = "row",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 4,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9,cutree_rows = 2)
h <- pheatmap(pd[which(apply(pd,1,max)>2),ord],scale = "row",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 5,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9,cutree_rows = 2)
h <- pheatmap(pd[which(apply(pd,1,max)>2),ord],scale = "row",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 4,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9,cutree_rows = 2)
h <- pheatmap(pd[which(apply(pd,1,max)>2),ord],scale = "row",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 4,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9,cutree_rows = 3)
rownames(pd.2)[h$tree_row$order]
which(rownames(pd.2)[h$tree_row$order]=='Podnl1')
sep <- which(rownames(pd.2)[h$tree_row$order]=='Podnl1')
h <- pheatmap(pd[which(apply(pd,1,max)>2),ord],scale = "row",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 4,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9,cutree_rows = 4)
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
p <- barfunMerge.2()
p
barfunMerge.2(sep.1 = sep+1,
              sep.2 = which(rownames(pd.2)[h$tree_row$order]=='Ccl1'))
barfunMerge.2(sep.1 = which(rownames(pd.2)[h$tree_row$order]=='Ccl1')+1,
              sep.2 = nrow(pd.2))
dim(NESvWT.glist)
ngeneFunc<- function(th.fdr,th.lfc) {
  edger.NESvWT.DEGs <- lapply(lrt,function(x)
    decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
  # NESvWT.glist
  g.cates <- calculate.overlap(lapply(edger.NESvWT.DEGs,
                                      function(x)which(x==-1)))
  NESvWT.glist <- data.frame(
    ids = unlist(g.cates),
    cates = unlist(sapply(1:7,function(x) rep(x,length(g.cates[[x]]))))
  )
  NESvWT.glist<- NESvWT.glist %>% mutate(ensemble = y$genes$Geneid[ids])
  NESvWT.glist<-NESvWT.glist%>% mutate(Symbol=getGsymble(ensemble))
  tmp.dic <- NESvWT.glist; rownames(tmp.dic) <- as.character(tmp.dic$ids)
  glist.NESred <- unique(unlist(lapply(edger.NESvWT.DEGs,
                                       function(x) which(x==-1))))
  glist.NESred.2 <- glist.NESred[glist.NESred%in% glist.wtup]
  pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
  sum(apply(pd,1,max)>=2)
}
ngeneFunc(1,0.05)
ngeneFunc(0.05,1)
ngene.fdr05 <- sapply(seq(0,1,by=.1),function(x) ngeneFunc(0.05,th.lfc = x))
ngeneFunc<- function(th.fdr,th.lfc) {
  edger.NESvWT.DEGs <- lapply(lrt,function(x)
    decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
  # NESvWT.glist
  g.cates <- calculate.overlap(lapply(edger.NESvWT.DEGs,
                                      function(x)which(x==-1)))
  NESvWT.glist <- data.frame(
    ids = unlist(g.cates),
    cates = unlist(sapply(1:7,function(x) rep(x,length(g.cates[[x]]))))
  )
  NESvWT.glist<- NESvWT.glist %>% mutate(ensemble = y$genes$Geneid[ids])
  #NESvWT.glist<-NESvWT.glist%>% mutate(Symbol=getGsymble(ensemble))
  tmp.dic <- NESvWT.glist; rownames(tmp.dic) <- as.character(tmp.dic$ids)
  glist.NESred <- unique(unlist(lapply(edger.NESvWT.DEGs,
                                       function(x) which(x==-1))))
  glist.NESred.2 <- glist.NESred[glist.NESred%in% glist.wtup]
  pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
  sum(apply(pd,1,max)>=2)
}
ngene.fdr05 <- sapply(seq(0,1,by=.1),function(x) ngeneFunc(0.05,th.lfc = x))
ngene.fdr05
pd.ngene.fdr05 <- data.frame(lfc= -seq(0,1,by=.1),
                             ngenes=sapply(seq(0,1,by=.1),function(x) ngeneFunc(0.05,th.lfc = x)))
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line()
?geom_text
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text()
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(lable=ngenes))
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes))
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =.5 )
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =1 )
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =2 )
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =4)
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =4,nudge_x = .05)
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =4,hjust=0)
ggplot(pd.ngene.fdr05,aes(lfc,ngenes)) + geom_point() +geom_line() +
  geom_text(aes(label=ngenes),nudge_y =4,hjust=1)
barfunMerge.2()

head(edger.NESvWT.DEGs$NESvWT_12_0)
colnames(edger.NESvWT.DEGs$NESvWT_12_0)
lrt$NESvWT_12_0$coefficients
head(lrt$NESvWT_12_0$coefficients)

## 2017-06-08
pheatmap(pd.lcpm[row_ord,],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = F,
         color= colorRampPalette(c("navy", "white", "firebrick1"))(10),
         breaks = seq(-m,m,length.out = 11),
         gaps_row = length(h1$tree_row$order),
         main="wt cat1&2 genes")
require(pheatmap)
pheatmap(pd.lcpm[row_ord,],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = F,
         color= colorRampPalette(c("navy", "white", "firebrick1"))(10),
         breaks = seq(-m,m,length.out = 11),
         gaps_row = length(h1$tree_row$order),
         main="wt cat1&2 genes")
load('./.RData')
ls()
names(edger.DEGs)
wt.genes <- list(A=which(edger.DEGs[[1]]==1&edger.DEGs[[2]]!=1),
                 B= which(edger.DEGs[[1]]==1&edger.DEGs[[2]]==1),
                 C= which(edger.DEGs[[1]]!=1&edger.DEGs[[2]]==1))
lapply(wt.genes,length)
pd <- lcpm[unlist(wt.genes),]
pd <- t(scale(t(pd))); m <- max(abs(range(pd)))
anno_row <- data.frame(category= sub("[0-9]+","",names(unlist(wt.genes))))
rownames(anno_row) <- rownames(pd)
ord <- (sample.info%>% mutate(ord=1:nrow(sample.info)) %>% arrange(desc(Genotype),
                                                                   treat))$ord
require(dplyr)
ord <- (sample.info%>% mutate(ord=1:nrow(sample.info)) %>% arrange(desc(Genotype),
                                                                   treat))$ord
h <- pheatmap(pd[anno_row=="A",ord],cluster_cols = F,scale = "none")
row_ord1 <- h$tree_row$order
h <- pheatmap(pd[anno_row=="B",ord],cluster_cols = F,scale = "none",
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              show_rownames = F,cutree_rows = 4)
row_ord2 <- h$tree_row$order + length(row_ord1)
h <- pheatmap(pd[anno_row=="C",ord],cluster_cols = F,scale = "none",
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21))
row_ord3 <- h$tree_row$order + length(row_ord2) + length(row_ord1)
row_ord <- c(row_ord1,row_ord2,row_ord3)
pheatmap(pd[row_ord,ord],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = F,
         color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
         breaks = seq(-m,m,length.out = 21),
         gaps_row = c(length(row_ord1),length(row_ord2)+length(row_ord1)),
         main="wt genes (only)",annotation_row = anno_row)
names(edger.DEGs)
glist.NESred.2 <- glist.NESred[glist.NESred%in% glist.wtup]
tmp.dic <- NESvWT.glist; rownames(tmp.dic) <- as.character(tmp.dic$ids)
pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
rownames(pd) <- tmp.dic[as.character(glist.NESred.2),]$Symbol
pd <- log2(pd+1)
pd <- t(scale(t(pd))); m <- max(abs(range(pd)))
ord <- c(1:3,7:9,13:15,4:6,10:12,16:18)
h <- pheatmap(pd[,ord],scale = "none",show_rownames = T,
              cluster_cols = F,cluster_rows = T,fontsize_row = 6,
              color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
              breaks = seq(-m,m,length.out = 21),border_color = NA,
              gaps_col = 9)
rownames(pd)[h$tree_row$order[1:13]]
length(glist.wtup)
length(glist.NESred)
length(glist.NESred.2)
th_fdr <- .05; th_lfc <- .5
ngeneFunc<- function(th.fdr,th.lfc) {
  edger.NESvWT.DEGs <- lapply(lrt,function(x)
    decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
  # NESvWT.glist
  g.cates <- calculate.overlap(lapply(edger.NESvWT.DEGs,
                                      function(x)which(x==-1)))
  NESvWT.glist <- data.frame(
    ids = unlist(g.cates),
    cates = unlist(sapply(1:7,function(x) rep(x,length(g.cates[[x]]))))
  )
  NESvWT.glist<- NESvWT.glist %>% mutate(ensemble = y$genes$Geneid[ids])
  #NESvWT.glist<-NESvWT.glist%>% mutate(Symbol=getGsymble(ensemble))
  tmp.dic <- NESvWT.glist; rownames(tmp.dic) <- as.character(tmp.dic$ids)
  glist.NESred <- unique(unlist(lapply(edger.NESvWT.DEGs,
                                       function(x) which(x==-1))))
  glist.NESred.2 <- glist.NESred[glist.NESred%in% glist.wtup]
  pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
  sum(apply(pd,1,max)>=2)
}
pd.ngene.fdr05 <- data.frame(lfc= -seq(0,1,by=.1),
                             ngenes=sapply(seq(0,1,by=.1),function(x) ngeneFunc(0.05,th.lfc = x)))
install.packages("Venndiagram")
install.packages("VennDiagram")
require(VennDiagram)
ngeneFunc<- function(th.fdr,th.lfc) {
  edger.NESvWT.DEGs <- lapply(lrt,function(x)
    decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
  # NESvWT.glist
  g.cates <- calculate.overlap(lapply(edger.NESvWT.DEGs,
                                      function(x)which(x==-1)))
  NESvWT.glist <- data.frame(
    ids = unlist(g.cates),
    cates = unlist(sapply(1:7,function(x) rep(x,length(g.cates[[x]]))))
  )
  NESvWT.glist<- NESvWT.glist %>% mutate(ensemble = y$genes$Geneid[ids])
  #NESvWT.glist<-NESvWT.glist%>% mutate(Symbol=getGsymble(ensemble))
  tmp.dic <- NESvWT.glist; rownames(tmp.dic) <- as.character(tmp.dic$ids)
  glist.NESred <- unique(unlist(lapply(edger.NESvWT.DEGs,
                                       function(x) which(x==-1))))
  glist.NESred.2 <- glist.NESred[glist.NESred%in% glist.wtup]
  pd <- rpkm[tmp.dic[as.character(glist.NESred.2),]$ensemble,]
  sum(apply(pd,1,max)>=2)
}
pd.ngene.fdr05 <- data.frame(lfc= -seq(0,1,by=.1),
                             ngenes=sapply(seq(0,1,by=.1),function(x) ngeneFunc(0.05,th.lfc = x)))
pd.ngene.fdr05 <- data.frame(lfc= -seq(0,1,by=.1),
                             ngenes=sapply(seq(0,1,by=.1),function(x) ngeneFunc(0.05,th.lfc = x)))
lfc= -seq(0,1,by=.1)
lfc
ls
ls()
wt.genes
dim(edger.DEGs)
dim(edger.DEGs$WT12_12)
length(unlist(wt.genes))
lapply(wt.genes,length)
wt.genes <- list(A=which(edger.DEGs[[1]]==1&edger.DEGs[[2]]!=1),
                 B= which(edger.DEGs[[1]]==1&edger.DEGs[[2]]==1),
                 C= which(edger.DEGs[[1]]!=1&edger.DEGs[[2]]==1))
lapply(wt.genes,length)
edger.DEGs
38+137+1372
require(venn)
install.packages()
install.packages('venn')
require(venn)
venn(NESvWT.glist.l2,zcolor = "style",cexil = 0.75)
dev.off()
dev.off()
venn(NESvWT.glist.l2,zcolor = "style",cexil = 0.75)
savehistory("~/Dropbox/Projects/MISC_bioinfo/T-Cell-project/Rhistory_2017-06-08.Rhistory")
