
# plot WT genes  ----------------------------------------------------------
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

# wt focus only WT --------------------------------------------------------
mt.genes.fil <- which(edger.DEGs[[3]]==1|edger.DEGs[[4]]==1)
wt.genes <- lapply(wt.genes,function(x) x[!x%in% mt.genes.fil])
lapply(wt.genes.2,length)

write.table(rownames(pd)[!rownames(pd)%in% rownames(lcpm)[unlist(wt.genes.2)]],
            file = 'wt_mt_overlap_gene.txt',quote = F,sep = "\n",row.names = F,col.names = F)
# save edgeR results 
pd <- lapply(1:2,function(i)
  data.frame(topTags(lrt[[i]],sort.by = 'none',n = nrow(lrt[[i]]))$table[,c(1,6:11)],
    compare=paste0(names(lrt)[i],".vs.WT12_0"),
    stringsAsFactors = F))
pd <- do.call(rbind,pd); 
pd$Symbol <- getGsymble(pd$Geneid)
write.csv(file='wt_edgeR_result.csv',pd,quote = F,row.names = F)
# wt LFC ------------------------------------------------------------------

pd <- lapply(1:4,function(i)
  topTags(lrt[[i]],sort.by = 'none',n = nrow(lrt[[i]]))$table[unlist(wt.genes),"logFC"])

pd <- do.call(cbind,pd); colnames(pd) <- names(lrt)
rownames(pd) <- rownames(lcpm[unlist(wt.genes),])

m <- 4; pd[pd>m] <- m ; pd[pd < -m] <- -m
#m <- ceiling(max(pd));
seps <- c(-rev(seq(1,m,by = .5)),0,seq(1,m,by = .5))
cols <- colorRampPalette(c("navy", "white", "firebrick3"))(length(seps)-1)
cols[c(length(cols)/2,length(cols)/2+1)] <- grey(.95)
pheatmap(pd[row_ord,],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = F,
         color=cols,
         
         breaks = seps,
         gaps_row = c(length(row_ord1),length(row_ord2)+length(row_ord1)),
         main="wt genes (all)",annotation_row = anno_row,border_color = NA)

# paried heatmap  ---------------------------------------------------------

require(LSD)
png(file='paried_scatter.png',width = 1600,height = 1600)
heatpairs(lcpm[,ord[h$tree_row$order]])
dev.off()


require(GGally)
png(file='paried_scatter2.png',width = 14,height = 14)
ggpairs(lcpm[,ord[h$tree_row$order]])
dev.off()

cor.mat <- cor(lcpm[,ord])
h<- pheatmap(cor.mat)

# mt only genes -----------------------------------------------------------
names(edger.DEGs)
mt.genes <- list(A=which(edger.DEGs[[3]]==1&edger.DEGs[[4]]!=1),
                 B= which(edger.DEGs[[3]]==1&edger.DEGs[[4]]==1),
                 C= which(edger.DEGs[[3]]!=1&edger.DEGs[[4]]==1))
lapply(mt.genes,length)
wt.genes.fil <- which(edger.DEGs[[1]]==1|edger.DEGs[[2]]==1)
mt.genes <- lapply(mt.genes,function(x) x[!x%in% wt.genes.fil])


pd <- lcpm[unlist(mt.genes),]
pd <- t(scale(t(pd))); m <- max(abs(range(pd)))
anno_row <- data.frame(category= sub("[0-9]+","",names(unlist(mt.genes))))
rownames(anno_row) <- rownames(pd)
ord <- (sample.info%>% mutate(ord=1:nrow(sample.info)) %>% arrange(desc(Genotype),
                                                                   treat))$ord

h <- pheatmap(pd[anno_row=="A",ord],cluster_cols = F,scale = "none")
row_ord1 <- h$tree_row$order
h <- pheatmap(pd[anno_row=="B",ord],cluster_cols = F,scale = "none")
row_ord2 <- h$tree_row$order + length(row_ord1)
h <- pheatmap(pd[anno_row=="C",ord],cluster_cols = F,scale = "none")
row_ord3 <- h$tree_row$order + length(row_ord2) + length(row_ord1)
row_ord <- c(row_ord1,row_ord2,row_ord3)


pheatmap(pd[row_ord,ord],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = F,
         color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
         breaks = seq(-m,m,length.out = 21),
         gaps_row = c(length(row_ord1),length(row_ord2)+length(row_ord1)),
         main="mt genes (only)",annotation_row = anno_row)

write.table(rownames(pd),
            file = 'kb_genes.txt',quote = F,sep = "\n",row.names = F,col.names = F)

# mt only genes lfc -------------------------------------------------------

pd <- lapply(1:4,function(i)
  topTags(lrt[[i]],sort.by = 'none',n = nrow(lrt[[i]]))$table[unlist(mt.genes),"logFC"])

pd <- do.call(cbind,pd); colnames(pd) <- names(lrt)
rownames(pd) <- rownames(lcpm[unlist(mt.genes),])


m <- ceiling(max(pd));
seps <- c(-rev(seq(1,m,by = .5)),0,seq(1,m,by = .5))
cols <- colorRampPalette(c("navy", "white", "firebrick3"))(length(seps)-1)
cols[c(length(cols)/2,length(cols)/2+1)] <- grey(.95)
pheatmap(pd[row_ord,],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = F,
         color=cols,
         
         breaks = seps,
         gaps_row = c(length(row_ord1),length(row_ord2)+length(row_ord1)),
         main="mt only genes ",annotation_row = anno_row,border_color = NA)

# wt - gene LFC vs. mt  ---------------------------------------------------
colnames(pd)
 pd.2 <-pd

pd <- as.data.frame(pd) %>% 
  mutate(ensemblID=rownames(pd.2))%>%
  rowwise%>%
  mutate(WT_mlfc = max(WT12_12,WT24_24),
         NES_mlfc = max(NES12_12,NES24_24))
pd <- pd %>% mutate(dmlfc=NES_mlfc-WT_mlfc)
pd$Symbol <- getGsymble(pd$ensemblID)

bks <- with(pd,c(min(dmlfc)-.0001,-1, -log2(1.5),log2(1.5),1,max(dmlfc)+0.001))
pd$cate=with(pd,cut(dmlfc,bks))
pd$cate2 <- as.numeric(pd$cate)

write.csv(file='wt_genes_lfc.csv',pd,quote = F,row.names = F)

table(pd$cate)
pd.3 <- pd; pd.3$cate<- as.numeric(pd$cate)

                                               
                                                                                           
#heatscatter(pd$WT_mlfc,pd$NES_mlfc,xlim = range(pd[,5:6]),
#            ylim=range(pd[,5:6]))
cols <- colorRampPalette(c("chartreuse4","white","brown"))(5)
plot(pd$WT_mlfc,pd$NES_mlfc,xlim = range(pd[,6:7]),
            ylim=range(pd[,6:7]),pch=21,bg=cols[pd$cate],
      xlab="Max WT LFC",ylab="Max NES LFC",main="kb genes")
sapply(bks[-c(1,length(bks))],function(x)
  abline(a=x,b=1,col=grey(.2),lty=2))
abline(a=0,b=1,col=grey(.2),lty=1)
grid()

lg <- table(pd$cate)
legend(x=0,y=10.3,pch = 21,pt.bg = cols,
       legend = paste(lg,"genes"),bg = n)




plot(pd$WT_mlfc,pd$NES_mlfc,xlim = range(pd[,5:6]),
     ylim=range(pd[,5:6]),pch=21,bg=cols[pd$cate2],
     xlab="Max WT LFC",ylab="Max NES LFC")
abline(a=-log2(1.5),b=1,col=grey(.2),lty=2)
abline(a=log2(1.5),b=1,col=grey(.2),lty=2)
abline(a=0,b=1,col=grey(.2),lty=1)
abline(v = 1,lty=2,col=grey(.8))
lg <- table(pd$cate2)
legend(x=0,y=10,pch = 21,col=cols,legend = paste(lg,"genes"))


# gage analysis -----------------------------------------------------------
require(gage)
kegg.mm <- kegg.gsets(species = "mouse")
exp.fc <- pd$NES_mlfc-pd$WT_mlfc;
nms <- names(exp.fc) <- rownames(pd)
require(biomaRt)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
nms.entrez.2 <- getBM(attributes=c('ensembl_gene_id','entrezgene'),
                    filters = 'ensembl_gene_id',
                    values  =nms, mart = mouse)


require(org.Mm.eg.db)
nms.entrez<- mapIds(org.Mm.eg.db,
       keys=nms,
       column="ENTREZID",
       keytype="ENSEMBL",
       multiVals="first")
names(exp.fc) <- nms.entrez

getPathId<-function(exp.fc){
  fc.kegg.p <- gage(2^exp.fc, gsets = kegg.mm$kg.sets,
                    ref = NULL, samp = NULL)
  sel <- fc.kegg.p$greater[, "q.val"] < 0.01 &
    !is.na(fc.kegg.p$greater[, "q.val"])
  path.ids <- rownames(fc.kegg.p$greater)[sel]
  sel.l <- fc.kegg.p$less[, "q.val"] < 0.01 &
    !is.na(fc.kegg.p$less[,"q.val"])
  path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
  path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)}


gage.test <- getPathId(exp.fc)

lcpm.2 <- lcpm 
rownames(lcpm.2) <-  mapIds(org.Mm.eg.db,
                            keys=rownames(lcpm),
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first")
ref.idx <- c(1,7,13); samp.idx <- c(5,11,17)
fc.kegg.p <- gage(lcpm.2,gsets = kegg.mm$kg.sets,ref = ref.idx,
                    samp = samp.idx,compare = "unpaired")

# plot 
substr(path.ids[2],1,8)

require(pathview)
cnts.d= lcpm.2[, samp.idx+1]-rowMeans(lcpm.2[, ref.idx+1])

pathview( gene.data = cnts.d, pathway.id = "mmu04141",
          gene.idtype = "KEGG",
    species = "mouse", same.layer = T)

write.csv(file='gage_kegg_up.csv', fc.kegg.p$greater[sel,],quote = F)
write.csv(file='gage_kegg_dwon.csv', fc.kegg.p$less[sel.l,],quote = F)


# Check kb genes ----------------------------------------------------------
nfkb.genes <- read.csv(file="nfkb.genes.txt",stringsAsFactors = F,
                       header = F)
nfkb.genes <- read.csv(file="nfkb_genes_kim.csv",stringsAsFactors = F,
                       header = T)
nfkb.genes<-nfkb.genes %>% filter(isCyto)%>% dplyr::select(Geneid,Symbol)


pd <- lcpm[rownames(lcpm)%in% nfkb.genes$V1,]
pd <- lcpm[rownames(lcpm)%in% nfkb.genes$Geneid,]
pd <- lapply(1:4,function(i)
  topTags(lrt[[i]],sort.by = 'none',
          n = nrow(lrt[[i]]))$table[rownames(lcpm)%in% nfkb.genes$Geneid,c("logFC")])
pd <- do.call(cbind,pd); 
colnames(pd) <- names(lrt)
rownames(pd) <- rownames(lcpm[rownames(lcpm)%in% nfkb.genes$Geneid,])
pd.2 <-pd
pd <- as.data.frame(pd) %>% 
  mutate(ensemblID=rownames(pd.2))%>%
  rowwise%>%
  mutate(WT_mlfc = max(WT12_12,WT24_24),
         NES_mlfc = max(NES12_12,NES24_24))
pd <- pd %>% mutate(dmlfc=NES_mlfc-WT_mlfc)
pd$Symbol <- getGsymble(pd$ensemblID)
pd$cate2 <- as.numeric(pd$cate)

plot(pd$WT_mlfc,pd$NES_mlfc,xlim=range(pd[,5:6 +1]),
     ylim=range(pd[,5:6+1]),
     pch=21,bg=cols[pd$cate2],
     xlab="WT max LFC",ylab="NES max LFC",
     main=paste(nrow(pd),"kb genes"))
lg <- table(pd$cate2)
legend(x=-1,y=5,pch = 21,pt.bg = cols,legend = paste(lg,"genes"))
abline(v=1,lty=2,col='grey');grid()
abline(a=0,b=1,col=grey(.2))
sapply(bks[-c(1,length(bks))],function(x)
  abline(a=x,b=1,col=grey(.2),lty=2,lwd=1))



bks <- with(pd,c(min(dmlfc)-.0001,-1, -0.5,0.5,1,max(dmlfc)+0.001))
pd$cate=with(pd,cut(dmlfc,bks))

table(pd$cate)
pd.3 <- pd; pd.3$cate<- as.numeric(pd$cate)


pheatmap(pd[,ord],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color= colorRampPalette(c("navy", "white", "firebrick3"))(20),
         breaks = seq(-m,m,length.out = 21),
         main="kb genes ")

write.csv(file='kb_genes_lfc.csv',pd,quote = F,row.names = F)
s.tmp <- pd$WT_mlfc;names(s.tmp)<- pd$ensemblID
barcodeplot(s.tmp,index =  pd$ensemblID%in% nfkb.genes$Geneid)
barcodeplot(pd.basal$NES_ctrl_LCPM - pd.basal$WT_ctrl_LCPM,
            index = rownames(pd.basal) %in% nfkb.genes$Geneid)
# check basal all  --------------------------------------------------------

pd.basal <- data.frame(WT_ctrl_LCPM=apply(lcpm[,grep("WT_12_0",colnames(lcpm))],1,mean),
                       NES_ctrl_LCPM=apply(lcpm[,grep("NES_12_0",colnames(lcpm))],1,mean))

pd.basal <- pd.basal[,-c(ncol(pd.basal)-1,ncol(pd.basal))]

heatscatter(pd.basal$WT_ctrl_LCPM,pd.basal$NES_ctrl_LCPM)
pd.basal.nfkb <- pd.basal[rownames(pd.basal)%in% nfkb.genes$Geneid,]
pd.basal.nfkb$cate2 <- as.numeric(pd.basal.nfkb$cate)
plot(pd.basal.nfkb$WT_ctrl_LCPM,pd.basal.nfkb$NES_ctrl_LCPM,
     xlim=range(pd.basal.nfkb[,1:2]),
     ylim=range(pd.basal.nfkb[,1:2]),
     pch=21,bg=cols[pd.basal.nfkb$cate2],
     xlab="WT basal LCPM",ylab="NES basal LCPM",
     main=paste(nrow(pd.basal.nfkb),"kb genes"))
lg <- table(pd.basal.nfkb$cate2)
legend(x=-4,y=10,pch = 21,pt.bg = cols,legend = paste(lg,"genes"))
grid()
abline(a=0,b=1,col=grey(.2))
sapply(bks[-c(1,length(bks))],function(x)
  abline(a=x,b=1,col=grey(.2),lty=2,lwd=1))

pd.basal.nfkb %>%filter(cate2==1) %>% arrange(WT_ctrl_LCPM)

with(pd.basal.nfkb,
     heatscatter(WT_ctrl_LCPM,NES_ctrl_LCPM,
                 xlim = range(c(WT_ctrl_LCPM,NES_ctrl_LCPM)),
                 ylim = range(c(WT_ctrl_LCPM,NES_ctrl_LCPM)),
                 main = paste(nrow(pd.basal.nfkb),"kb genes")))
abline(a=0,b=1,col=grey(0.2))
sapply(bks[-c(1,length(bks))],function(x)
  abline(a=x,b=1,col=grey(.2),lty=2,lwd=1))


pd.basal.nfkb$ensembl <- rownames(pd.basal.nfkb)
pd.basal.nfkb$symbol <- getGsymble(rownames(pd.basal.nfkb))
pd.basal.nfkb$delta_lcpm <- pd.basal.nfkb$NES_ctrl_LCPM -pd.basal.nfkb$WT_ctrl_LCPM
pd.basal.nfkb$cate <- with(pd.basal.nfkb,
                           cut(delta_lcpm,c(min(delta_lcpm)-0.001,-1,-0.5,0.5,1,max(delta_lcpm)+0.001)))
table(pd.basal.nfkb$cate)

pd.basal.nfkb%>%
  filter(cate=="(-2.53,-1]"|cate=="(-1,-0.5]") %>%
  arrange(WT_ctrl_LCPM)

heatscatter(lcpm[,1],lcpm[,7])
t.test(x=pd.basal[,1],y=pd.basal[,2],paired = T,alternative = "greater")
t.test(x=lcpm[,1],y=lcpm[,7],paired = T,alternative = "greater")




heatscatter(x=pd.basal[unlist(wt.genes),1],pd.basal[unlist(wt.genes),2])
abline(lm(pd.basal[unlist(wt.genes),2]~pd.basal[unlist(wt.genes),1]),col=2)


# hist 
pd.basal.wt <- pd.basal[unlist(wt.genes),]
pd.basal.wt$ensembl <- rownames(pd.basal.wt)
pd.basal.wt$symbol <- getGsymble(pd.basal.wt$ensembl)
pd.basal.wt$delta <- pd.basal.wt$NES_ctrl_LCPM-pd.basal.wt$WT_ctrl_LCPM
pd <- pd.basal.wt
pd$cate <- with(pd,cut(delta,c(min(delta)-0.001,-1,-0.5,0.5,1,max(delta)+0.001)))
pd$cate2 <- as.numeric(pd$cate)
pd.basal.wt$cate <- pd$cate
p.1 <- ggplot(pd.basal.wt,aes(WT_ctrl_LCPM)) +
  geom_density(aes(colour=cate)) + scale_color_manual(values = cols)
p.2 <- ggplot(pd.basal.wt,aes(NES_ctrl_LCPM)) + 
  geom_density(aes(colour=cate))+
  scale_color_manual(values = cols)

pd.2 <- read.csv(file="wt_genes_lfc.csv",stringsAsFactors = F,
                 head = F,skip = 1)
colnames(pd.2) <- c(read.csv(file="wt_genes_lfc.csv",stringsAsFactors = F,
                           head = F,nrows = 1),"cate2.2")
pd.3 <- pd
pd <- pd.3[pd.2$cate2.2%in% c(1,2),]
plot(pd$WT_ctrl_LCPM,pd$NES_ctrl_LCPM,xlim=range(pd[,1:2]),
     ylim=range(pd[,1:2]),
     pch=21,bg=cols[pd.2$cate2.2[pd.2$cate2.2%in% c(1,2)]],
     xlab="WT basal LCPM",ylab="NES basal LCPM",
     main=paste(nrow(pd)," genes"))
lg <- table(pd$cate2)
legend(x=-5,y=6,pch = 21,pt.bg = cols,legend = paste(lg,"genes"))
grid()
abline(a=0,b=1,col=grey(.2))
sapply(bks[-c(1,length(bks))],function(x)
  abline(a=x,b=1,col=grey(.2),lty=2,lwd=1))


head(lcpm[pd.3$ensemblID[pd.3$cate==1],])

write.csv(file='kb_genes_basal.csv',pd.basal.nfkb,quote = F,row.names = F)

# pheatmap 

pd.lcpm <- rbind(lcpm[pd[pd$cate2==1,]$ensemblID,ord],
                 lcpm[pd[pd$cate2==2,]$ensemblID,ord])
pd.lcpm <- t(scale(t(pd.lcpm)))
m<-max(abs(pd.lcpm))
h1 <- pheatmap(lcpm[pd[pd$cate2==1,]$ensemblID,ord],scale = 'row',
               cluster_cols = F)
h2 <- pheatmap(lcpm[pd[pd$cate2==2,]$ensemblID,ord],scale = 'row',
               cluster_cols = F)
row_ord1 <- h1$tree_row$order
row_ord2 <- h2$tree_row$order+length(row_ord1)

pheatmap(pd.lcpm[row_ord,],scale = "none",show_rownames = F,
         cluster_cols = F,cluster_rows = F,
         color= colorRampPalette(c("navy", "white", "firebrick1"))(10),
         breaks = seq(-m,m,length.out = 11),
         gaps_row = length(h1$tree_row$order),
         main="wt cat1&2 genes")
