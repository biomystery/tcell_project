require(edgeR)
require(tidyverse)
require(ggplot2)
require(pheatmap)

pd <- list()
geom_noboarder <- theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))
getGsymble <- function(glist){
  require(biomaRt)
  require(org.Mm.eg.db)
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

plotLegend <- function(cols,bks,fnames){
  require(RColorBrewer)
  #bks <- seq(round(min(rsums)),round(max(rsums))+1,length.out = 6)
  #cols<- colorRampPalette(c( "white", "blueviolet"))(5)
  #fnames<-'tmp.eps'
  setEPS()
  postscript(fnames,onefile = F,width = 0.1,height = .05*2*length(bks))
  par(mar=c(0.02, 0.04, 0.04,0.1))
  barplot(rep(1,length(bks)-1),width=diff(bks),space = 0,border = NA,
          col = cols,axes = F,horiz = T,
          xaxs = "i", yaxs = "i",xlim = c(0,1),
          ylim=c(0,sum(diff(bks)))
          #xlim=c(0,1)
  )
  #axis(4,labels = F,tck = 0.2,at=cumsum(diff(bks)))
  #box()
  dev.off()
}
# Fig.5A -maPlot ----------------------------------------------------------

### data generate
if(T){
  ## read counts & sample info
  raw.count <- read.delim(file="./data/counts-gene.txt",header = T,skip = 1,
                          stringsAsFactors = F)
  anno.count <- raw.count[,1:6]; raw.count <- raw.count[,-c(1:6)]
  colnames(raw.count) <- sub(".filtered.bam","",colnames(raw.count)) 
  sample.info <- read.csv(file="./data/sample_info.csv",stringsAsFactors = F,
                          header = T)
  anno.count$Geneid <- substr(anno.count$Geneid,1,18)
  sample.info$Number <- colnames(raw.count)
  sample.info$treat<-with(sample.info,paste(Genotype,X_CD3.CD28..hrs.,X_4.1BB..hrs.,sep = "_"))
  sample.info$treat <- factor(sample.info$treat)
  
  ## construct DGElist object 
  y <- DGEList(counts = raw.count,genes = anno.count)
  y <- y[rowSums(cpm(y)>1)>=2,,keep.lib.sizes=F] #at least two sample cpm>2,11664x18
  y <- calcNormFactors(y)
  sample.info$treat<-with(sample.info,paste(X_CD3.CD28..hrs.,X_4.1BB..hrs.,sep = "_"))
  design <- model.matrix(~0+ Genotype:treat ,sample.info)
  y <- estimateDisp(y,design)
  
  ## glmFit 
  fit <- glmFit(y,design )
  
  ## compare vs. t=0
  edger.comp.mat <- data.frame(WT12_12= c(0,-1,0,1,0,0),
                               WT24_24=c(0,-1,0,0,0,1),
                               NES12_12=c(-1,0,1,0,0,0),
                               NES24_24=c(-1,0,0,0,1,0))
  lrt <- lapply(colnames(edger.comp.mat),function(x)
    glmLRT(fit,contrast = edger.comp.mat[,x]))
  th.fdr <- .01; th.lfc <- 1
  edger.DEGs <- lapply(lrt,function(x) decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
  
  ## compare nes vs. WT
  edger.comp.mat <- data.frame(
    NESvWT_12_0= c(1,-1,0,0,0,0),
    NESvWT_12_12= c(0,0,1,-1,0,0),
    NESvWT_24_24=c(0,0,0,0,1,-1))
  
  lrt.NESvWT <- sapply(colnames(edger.comp.mat),function(x)
    glmLRT(fit,contrast = edger.comp.mat[,x]),simplify = F)
  
  th.fdr <- .05; th.lfc <- .5
  edger.NESvWT.DEGs <- lapply(lrt.NESvWT,function(x) decideTestsDGE(x,p.value = th.fdr,lfc = th.lfc))
}

## subfig5a_wt
if(T){
  setEPS()
  postscript(file="./subfig5a_wt.eps",width = 1.25,height = 1)
  par(mar=c(0.1,.1,.1,.1))
  i=1
  plotSmear(lrt[[i]],de.tags = rownames(y)[as.logical(edger.DEGs[[i]])],
            deCol="#ed472e",smooth.scatter = T)
  abline(h=c(-1, 1), col="blue")
  #text(10,4,paste0('Up:',summary(edger.DEGs[[i]])[3],'\n Down:',summary(edger.DEGs[[i]])[1]))
  
  dev.off()
}

## subfig5a_nes
if(T){
  setEPS()
  postscript(file="./subfig5a_nes.eps",width = 1.25,height = 1)
  par(mar=c(0.1,.1,.1,.1))
  i=2
  #plotSmear(lrt.NESvWT[[i]],de.tags = rownames(y)[as.logical(edger.NESvWT.DEGs[[i]])],
  #          main=colnames(edger.comp.mat)[i])
  plotSmear(lrt.NESvWT[[i]],de.tags = rownames(y)[as.logical(edger.NESvWT.DEGs[[i]])],
            deCol='#ed472e',smooth.scatter = T)
  
  abline(h=c(-th.lfc, th.lfc), col="blue")
  #text(10,3,paste0('Up:',summary(edger.NESvWT.DEGs[[i]])[3],'\n Down:',summary(edger.NESvWT.DEGs[[i]])[1]))
  dev.off()
}

### save data 
pd$a_mdplot$y =y;pd$a_mdplot$fit =fit;
pd$a_mdplot$lrt =lrt;  pd$a_mdplot$edger.DEGs = edger.DEGs
pd$a_mdplot$lrt.NESvWT =lrt.NESvWT;  pd$a_mdplot$edger.DEGs = edger.DEGs

### - venn diagam 
if(T){
  pd$b_venn$id_list=list(wt_induced = which(edger.DEGs[[1]]==1),
                         nes_reduced =which(edger.NESvWT.DEGs[[2]]==-1))
  pd$b_venn$nes_idx = intersect(pd$b_venn$id_list$wt_induced,pd$b_venn$id_list$nes_reduced)
  pd$a_mdplot$venn$nes_ind_idx <- setdiff(pd$a_mdplot$venn$id_list$wt_induced,pd$a_mdplot$venn$id_list$nes_reduced)
  pd$a_mdplot$venn$nes_red_idx <- setdiff(pd$a_mdplot$venn$id_list$nes_reduced,pd$a_mdplot$venn$id_list$wt_induced)
  require(venn)
  venn(pd$b_venn$id_list,
       zcolor = "style",ellipse = F,cexil = 0.75)
}

# Fig. 5B kegg_pathway analysis  ----------------------------------------------------------

### generate bg and save bg and nes glist

#set.seed(111)
#bg_list <- sample(y$genes$Geneid[-pd$b_venn$nes_idx],length(pd$b_venn$nes_idx))
#write.table(file='./finalize/KEGG_pathway/bg_list.txt',bg_list,quote = F,row.names = F,col.names = F)
write.table(file='./finalize/KEGG_pathway/nes_list.txt',y$genes$Geneid[pd$b_venn$nes_idx],
          quote = F,row.names = F,col.names = F)
pd.kegg <- read.delim(file = './finalize/KEGG_pathway/david_go_nes_03_31_2018.txt')
pd.kegg <- pd.kegg%>%filter(PValue<0.05) %>% arrange(Fold.Enrichment)
pd.kegg$Term<-sub("mmu[0-9]+:","",pd.kegg$Term)
pd.kegg$Term <- factor(pd.kegg$Term,levels = pd.kegg$Term)

### plot
p <- ggplot(pd.kegg,aes(y=Fold.Enrichment,x=Term))+ geom_bar(stat = 'identity',fill="#ed472e")+ coord_flip()+
  theme_bw()+theme(axis.title = element_blank(),text = element_text(size = 15))
ggsave(filename = "Subfig7c.pdf",width = 3,height = 2.5,plot = p,
       scale = 2.5)
ggsave(filename = "Subfig7c.eps",width = 1.2,height = 2.4,plot = p+theme(axis.text = element_blank()),
       scale = 2.5)
pd$a_mdplot$venn <- pd$b_venn
pd$b_venn <- NULL
pd$b_kegg <- pd.kegg

# Fig.5c heatmap  ---------------------------------------------------------

## generat rpkm & data 
dat.rpkm <- rpkm(y)
dat.hm <- lapply(list(nes=pd$a_mdplot$venn$nes_idx,
                      nes_ind=pd$a_mdplot$venn$nes_ind_idx,
                      nes_red=pd$a_mdplot$venn$nes_red_idx
                      ),function(x){
  tmp <- dat.rpkm[x,c(1,7,13,3,9,15,5,11,17,2,8,14,4,10,16,6,12,18)]
  tmp <- tmp[,-grep("24_24",colnames(tmp))]
  nm <- getGsymble(y$genes$Geneid[x])
  rownames(tmp) <- as.character(nm)
  tmp
})

ll=list(nes=pd$a_mdplot$venn$nes_idx,
     nes_ind=pd$a_mdplot$venn$nes_ind_idx,
     nes_red=pd$a_mdplot$venn$nes_red_idx)
for(i in 1:3){
  write.csv(cbind(dat.hm[[i]],ensembl_id=y$genes$Geneid[ll[[i]]]),
            file = paste0("Table_",names(dat.hm)[i],"_rpkm.csv"),quote = F)
}                      
          

### ord hm - get order
tmp <- pd$c_hm
pd$c_hm <- pd$c_hm[c("nes_ind","nes","nes_red")]

id.list <- lapply(pd$c_hm, function(x){
  pheatmap(x,scale = "row",
           cluster_cols = F,cluster_rows = T)
})

dat.hm <- lapply(names(pd$c_hm), function(nm){
  pd$c_hm[[nm]][id.list[[nm]]$tree_row$order,]
})
dat.hm <- do.call(rbind,dat.hm)
mn <- -max(abs(range(t(scale(t(dat.hm))))));mx<- -mn
dat.hm.2 <- t(scale(t(dat.hm)))
dat.hm.2[dat.hm.2>2]<- 2; dat.hm.2[dat.hm.2< -2]<- -2;m = 2.01;

# or 0-1 
dat.hm.2 <- t(apply(dat.hm,1, norm.func))
mn<- 0; mx<-1
if(T){
  n=50;
  setEPS()
  postscript(file="./subfig5c_hm_anno.eps",width = 5,height = 12)
  pheatmap(dat.hm.2,scale = "none",cluster_cols = F,cluster_rows = F,
           gaps_row = cumsum(sapply(pd$c_hm,nrow)),fontsize_row = 3,
           breaks = seq(mn,mx,length.out = n+1),show_colnames =F,
           color= colorRampPalette(c("white", "#ed472d"))(n),
           gaps_col = c(3,6,9))
  dev.off()
  setEPS()
  postscript(file="./subfig5c_hm.eps",width = 2.4,height = 3)
  pheatmap(dat.hm.2,scale = "none",cluster_cols = F,cluster_rows = F,
           gaps_row = cumsum(sapply(pd$c_hm,nrow)),show_rownames = F,
           breaks = seq(-mn,mx,length.out = n+1),show_colnames =F,
           color= colorRampPalette(c("white", "#ed472d"))(n), #firebrick3,navy
           #color= colorRampPalette(brewer.pal(9,"Reds"))(n),
           gaps_col = c(3,6,9),legend = F)
  dev.off()
 
  plotLegend(cols =colorRampPalette(c("white", "#ed472d"))(n),
            bks=seq(-mn,mx,length.out = n+1),
            fnames = "./subfig5c_hm_scale.eps")
}



#dat.hm <- t(apply(dat.hm, 1,norm.func))
#m <- max(abs(range(t(scale(t(dat.hm))))))
pheatmap(dat.hm,scale = "row",treeheight_row = 0,show_rownames =T,show_colnames = F,
         cluster_cols = F,cluster_rows = T,fontsize_row = 4,#fontsize_col = 4,
         color= colorRampPalette(c("navy","white", "firebrick3"))(20),
         cellwidth = 5,cellheight = 4,
         breaks = seq(-m,m,length.out = 21),border_color = NA,#cellwidth = 2,
         gaps_col = c(3,6,9),legend = F)
dev.off()


### avg over reps 
pd$c_hm_avg<- lapply(pd$c_hm, function(x){
  tmp<-t(do.call(rbind,lapply(c(1,4,7,10),function(i) apply(x[,i:(i+2)],1,mean))))
  names(tmp)<- c("WT_0hr","WT_12hr","NES_0hr","NES_12hr")
  tmp
})
sapply(names(pd$c_hm_avg), function(x){
  dat.hm <- pd$c_hm_avg[[x]]
  #m <- max(abs(range(t(scale(t(dat.hm))))))
  dat.hm <- apply(dat.hm, 1,norm.func)
  ## no anno
  setEPS()
  postscript(file=paste0("./subfig5c_hm_",x,"_anno.eps"),width = 8,height = 5)
  pheatmap(dat.hm,scale = "none",
           cluster_cols = T,cluster_rows = F,fontsize_col = 6,#fontsize_col = 4,
           color= colorRampPalette(c( "white", "firebrick3"))(20),
           width = 3.125,height = 2.6597,
           breaks = seq(0,1,length.out = 21),border_color = NA)
  dev.off()
  ## no anno
  
  setEPS()
  postscript(file=paste0("./subfig5c_hm_",x,".eps"),width = 6,height = 2)
  pheatmap(dat.hm,scale = "none",treeheight_col = 0,show_rownames =F,show_colnames = T,
           cluster_cols = T,cluster_rows = F,fontsize_col = 4,#fontsize_col = 4,
           color= colorRampPalette(c("white", "firebrick3"))(20),
           cellwidth = 5,cellheight = 4,
           breaks = seq(0,1,length.out = 21),border_color = NA,#cellwidth = 2,
           legend = F)
  dev.off()
})


### 
pd$c_hm_avg_lfc = log2(pd$c_hm_avg+1)

# Fig.5EF-motifs ----------------------------------------------------------

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



# Fig.5EF-motifs-bar ----------------------------------------------------------

dat <- read.csv(file = "./analysis3/wholeMM10/all.anno.reduced.csv",stringsAsFactors = T)
wt_remain <- y$genes$Geneid[setdiff(pd$a_mdplot$venn$id_list$wt_induced,pd$a_mdplot$venn$nes_idx)]
nes_dep <- y$genes$Geneid[intersect(pd$a_mdplot$venn$id_list$wt_induced,pd$a_mdplot$venn$nes_idx)]

dat <- dat%>%mutate(gcate = ifelse(Nearest.Ensembl %in% nes_dep,"NES-dep",
                                   ifelse(Nearest.Ensembl %in% wt_remain, "NES-indep",
                                   "Background")))%>%
  distinct()

## rm duplidated
dat <- dat %>% filter(Nearest.Ensembl!="")
dat$Nearest.Ensembl <- as.character(dat$Nearest.Ensembl)
dup.id <- dat$Nearest.Ensembl[duplicated(dat$Nearest.Ensembl)]
dat.dup <- dat%>% filter(Nearest.Ensembl %in% dup.id) %>% arrange(Nearest.Ensembl)
dat.dup.red <- dat.dup[1,]
for(g in dup.id){
  idx <- which(dat.dup$Nearest.Ensembl==g)
  new.row <-dat.dup[idx[1],]
  for (i in 4:6) new.row[i] <- max(dat.dup[idx,i])
  for (i in 8:10) new.row[i] <- (sum(dat.dup[idx,i])>0)
  dat.dup.red <- rbind(dat.dup.red,new.row) 
}
dat.dup.red <- dat.dup.red[-1,]
dat <- dat%>% filter(!Nearest.Ensembl %in% dup.id)
dat <- rbind(dat,dat.dup.red)

dat <- dat[!duplicated(dat),]
rownames(dat)<- dat$Nearest.Ensembl

if(T){
  #th <-c(9.68472745523327,8.30143760313367,5.53485789893449)
  th <- rep(7,3)
  dat <- dat %>%   mutate(p65_p50=p65_p50_mouse.Best.Motif.log.odds.Score>=th[1],
                          crel_p50=crel_p50_mouse.Best.Motif.log.odds.Score>=th[2],
                          crel=crel_mouse.Best.Motif.log.odds.Score>=th[3])
  
  dat.sum <- dat %>% group_by(gcate)%>%
    summarise(crel_p50 = sum(crel_p50)/n()*100,
              p65_p50 = sum(p65_p50)/n()*100,
              crel = sum(crel)/n()*100)%>% 
    gather(key = "TF",value = "percent",c(2:4))
  
  dat.sum$gcate <- factor(dat.sum$gcate,levels = c( "Background", "NES-indep","NES-dep"))
  cols <- c("grey",rep('black',2))
  ggplot(dat.sum, aes(x=gcate,y=percent,fill=gcate))+ geom_bar(stat = 'identity',width =.8)+
    scale_fill_manual(values = cols) + facet_wrap(~ TF) + ylim(0,53)
  ggsave(filename = "./subfig5d_anno.pdf")
  sapply(unique(dat.sum$TF),function(x){
    p <- ggplot(dat.sum%>% filter(TF==x), aes(x=gcate,y=percent,fill=gcate))+ geom_bar(stat = 'identity',width =.8)+
      scale_fill_manual(values = cols) + ylim(0,53) +
     theme_bw() + theme(title = element_blank(),legend.position = "none",
                                            axis.text =  element_blank())
    ggsave(filename = paste0("./subfig4s_motif_",x,".eps"),p+geom_noboarder,
           width = 1,height = 1.2,scale = 2.5)
  })  
  
  ## test 
  test.res <- list()
  calFE = function(tb) tb[1,1]/(sum(tb[1,1],tb[2,1]))*sum(tb[1,2]+tb[2,2])/tb[1,2]
  ## crel
  test.res$crel<-sapply(c("NES-indep","NES-dep"),function(x){
    tb <- table(dat%>% filter(gcate %in% c(x,"Background"))%>%
                  dplyr::select(crel,gcate))
    c(fisher.test(tb[c("TRUE","FALSE"),c(x,"Background")],alternative="greater"),
      FE = calFE(tb[c("TRUE","FALSE"),c(x,"Background")]))
  })
  
 
  ## crel_p50
  test.res$crel_p50<-sapply(c("NES-indep","NES-dep"),function(x){
    tb <- table(dat%>% filter(gcate %in% c(x,"Background"))%>%
                  dplyr::select(crel_p50,gcate))
    c(fisher.test(tb[c("TRUE","FALSE"),c(x,"Background")],alternative="greater"),
      FE = calFE(tb[c("TRUE","FALSE"),c(x,"Background")]))
  })
  
  ## p65_p50
  test.res$p65_p50<-sapply(c("NES-indep","NES-dep"),function(x){
    tb <- table(dat%>% filter(gcate %in% c(x,"Background"))%>%
                  dplyr::select(p65_p50,gcate))
    c(fisher.test(tb[c("TRUE","FALSE"),c(x,"Background")],alternative="greater"),
      FE = calFE(tb[c("TRUE","FALSE"),c(x,"Background")]))
  })

  test.res
}



if(T){
  #test.res <- list()
  calFE = function(tb) tb[1,1]/(sum(tb[1,1],tb[2,1]))*sum(tb[1,2]+tb[2,2])/tb[1,2]
  calLB = function(tb) exp(log(tb[1,1]*tb[2,2]/tb[1,2]/tb[2,1])-
                             1.96*sqrt(sum(1/tb)))
  calUB = function(tb) exp(log(tb[1,1]*tb[2,2]/tb[1,2]/tb[2,1])+
                             1.96*sqrt(sum(1/tb)))
  
  ## crel

  test.res$crel.2<-sapply(c("NES-dep"),function(x){
    tb <- table(dat%>% filter(gcate %in% c(x,"NES-indep"))%>%
                  dplyr::select(crel,gcate))
    c(fisher.test(tb[c("TRUE","FALSE"),c(x,"NES-indep")],alternative="greater"),
      FE = calFE(tb[c("TRUE","FALSE"),c(x,"NES-indep")]))
  })
  
  
  test.res$crel_p50.2<-sapply(c("NES-dep"),function(x){
    tb <- table(dat%>% filter(gcate %in% c(x,"NES-indep"))%>%
                  dplyr::select(crel_p50,gcate))
    c(fisher.test(tb[c("TRUE","FALSE"),c(x,"NES-indep")],alternative="greater"),
          FE = calFE(tb[c("TRUE","FALSE"),c(x,"NES-indep")]))
  })
  
  
  test.res$p65_p50.2<-sapply(c("NES-dep"),function(x){
    tb <- table(dat%>% filter(gcate %in% c(x,"NES-indep"))%>%
                  dplyr::select(p65_p50,gcate))
    c(fisher.test(tb[c("TRUE","FALSE"),c(x,"NES-indep")],alternative="greater"),
          FE = calFE(tb[c("TRUE","FALSE"),c(x,"NES-indep")]))
  })
  
  cols <- c(p65_p50="#d8e6a8",cRel_p50="#c9e6ce",cRel_cRel="#9fdded")
  dat.motif <- data.frame(motif=c("p65_p50","cRel_p50","cRel_cRel"),
                          odds_ratio=c(test.res$p65_p50.2[[3]],
                                       test.res$crel_p50.2[[3]],
                                       test.res$crel.2[[3]]),
                          FE = c(test.res$p65_p50.2[[8]],
                                 test.res$crel_p50.2[[8]],
                                 test.res$crel.2[[8]]),
                          pval=signif(c(test.res$p65_p50.2[[1]],
                                        test.res$crel_p50.2[[1]],
                                        test.res$crel.2[[1]]),3))
  
  p<- ggplot(dat.motif,aes(motif,odds_ratio,fill=motif))+
    geom_bar(stat='identity') + scale_fill_manual(values = cols)+#coord_cartesian(ylim=c(1.6,2.2))+
    coord_flip(ylim=c(1.6,2.4))+theme_bw() 
  p +geom_text(aes(label=pval))
  
  ggsave(p +geom_text(aes(label=pval)),filename = "./subfig5_motif_all_anno.eps")
  ggsave(p+theme(text = element_blank(),legend.position = "none")+geom_noboarder,
         filename = "./subfig5_motif_all.eps",
         width = 2.4,height = 2)
  
}




## save data 


ll=list(nes=pd$a_mdplot$venn$nes_idx,
        nes_ind=pd$a_mdplot$venn$nes_ind_idx,
        nes_red=pd$a_mdplot$venn$nes_red_idx)
for(i in 1:3){
  write.csv(cbind(dat.hm[[i]],ensembl_id=y$genes$Geneid[ll[[i]]],
                  dat[y$genes$Geneid[ll[[i]]],c(4:6,8:10)]),
            file = paste0("Table_",names(dat.hm)[i],"_rpkm.csv"),quote = F)
}                      





#p.2<- ggplot(dat.motif,aes(motif,FE,fill=motif))+
#  geom_bar(stat='identity') + scale_fill_manual(values = cols)+
#  coord_flip()+theme_bw()

# Fig4s_a_go_bp -----------------------------------------------------------
pd.go <- read.delim(file = './finalize/KEGG_pathway/david_go_BP_nes_04_01_2018.txt')
pd.go.back <- pd.go
pd.go <- pd.go.back%>%filter(Bonferroni<0.05) %>% arrange(Fold.Enrichment)
pd.go$Term<-sub("GO:[0-9]+~","",pd.go$Term)
pd.go$Term <- factor(pd.go$Term,levels = pd.go$Term)

### plot
p<- ggplot(pd.go,aes(y=Fold.Enrichment,x=Term))+ geom_bar(stat = 'identity',fill="#ed472e")+ coord_flip()+
  theme_bw()+theme(axis.title = element_blank(),text = element_text(size = 15))
ggsave(filename = "Subfigs4A_go_anno.pdf",width = 5,height = 2.5,plot = p,
       scale = 2.5)
ggsave(filename = "Subfig7c.eps",width = 1.2,height = 2.4,plot = p+theme(axis.text = element_blank()),
       scale = 2.5)
pd$a_mdplot$venn <- pd$b_venn
pd$b_venn <- NULL
pd$b_kegg <- pd.go


# Fig4s_reps_scatter ------------------------------------------------------

pd.log2 <- log2(dat.rpkm+1)

if(T){
  #pdf(file = "fig4s_reps_scatter.pdf",height = 8,width = 6)
  #setEPS()
  pdf(file = "fig4s_reps_scatter.pdf",height = 8,width = 6)
  
  par(mfrow=c(4,3))
  sapply(c(2,4,1,3),function(m){
    idx <- as.numeric(which(design[,m]==1))
    apply(combn(idx,2),2,function(x) {
      i=x[1];j=x[2]
      smoothScatter(pd.log2[,x[1]],pd.log2[,x[2]],
                    xlab = colnames(pd.log2)[i],
                    ylab=colnames(pd.log2)[j],
                    colramp = colorRampPalette(c("white","#000099", "#00FEFF", "#45FE4F", 
                                                 "#FCFF00", "#FF9400", "#FF3100")),
                    main=paste("Pearson's ",expression(rho),"=",signif(cor(pd.log2[,x[1]],
                                                                           pd.log2[,x[2]]),3)))
      abline(a=0,b=1)  
    })
    
  })
  
  par(mfrow=c(1,1))
  dev.off()
}

## linear scale
pd.log2 <- dat.rpkm
if(T){
  pdf(file = "fig4s_reps_scatter.pdf",height = 8,width = 6)
  par(mfrow=c(4,3))
  sapply(c(2,4,1,3),function(m){
    idx <- as.numeric(which(design[,m]==1))
    apply(combn(idx,2),2,function(x) {
      i=x[1];j=x[2]
      plot(pd.log2[,i],pd.log2[,j],
                    xlab = colnames(pd.log2)[i],
                    ylab=colnames(pd.log2)[j],
                    cex=0.2,xlim=c(0,100),ylim=c(0,100),
                    main=paste("Pearson's ",expression(rho),"=",signif(cor(pd.log2[,x[1]],
                                                                           pd.log2[,x[2]]),3)))
      abline(a=0,b=1)  
    })
    
  })
  
  par(mfrow=c(1,1))
  dev.off()
}


