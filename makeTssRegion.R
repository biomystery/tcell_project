genes.anno <- read.delim(file='~/Dropbox/Projects/RelA_project/data/2016-01-29-RNASeq-Rec-b2/counts.txt',stringsAsFactors = F)
rownames(genes.anno) <- genes.anno$Geneid; genes.anno$Geneid <- NULL 

getTS <- function(anno,type='TSS'){
  s <- anno$Strand
  mutlipleS<- ifelse(type=='TSS',anno$Start,anno$End) 
  tmp <- unlist(strsplit(mutlipleS,split = ';'))
  ifelse(s=='+',tmp[1],tmp[length(tmp)])
}


tss.genes.anno <- data.frame(Strand=unlist(lapply(genes.anno$Strand,
                                                  function(x) unlist(strsplit(x,split = ';'))[1])))
#rownames(tss.genes.anno)<-rownames(genes.anno)
tss.genes.anno$TSS <- sapply(1:nrow(genes.anno),
                              function(x) as.numeric(getTS(genes.anno[x,],type='TSS')))
tss.genes.anno$Chr <- unlist(lapply(genes.anno$Chr,
                                     function(x) unlist(strsplit(x,split = ";"))[1]))

tss.genes.bed <- data.frame(Chr=tss.genes.anno$Chr,
                            Start=ifelse(tss.genes.anno$Strand=="+",tss.genes.anno$TSS-1000,
                                   tss.genes.anno$TSS-500),
                            End=ifelse(tss.genes.anno$Strand=="+",tss.genes.anno$TSS+500,
                                         tss.genes.anno$TSS+1000),
                            Strand = tss.genes.anno$Strand,
                            Ensembl=substr(rownames(genes.anno),1,18))

write.table(tss.genes.bed,file='mm10.promoter.bed',quote = F,sep = "\t",
            row.names = F,col.names = F)
