## ------------------------------------------------------------------------
library(rtracklayer)
library(ggplot2)
library(plyr)
library(reshape2)
library(plotly)
library(VennDiagram)
require(gplots) 
system.time(gtf_gencode_comp <- readGFF("~/Documents/HCA/reference/refs/gencode.v27.chr_patch_hapl_scaff.annotation.gtf", version=2L, tags = c("gene_name","gene_id", "transcript_id","gene_type")))
system.time(gtf_gencode_basic <- readGFF("~/Documents/HCA/reference/refs/gencode.v27.chr_patch_hapl_scaff.basic.annotation.gtf", version=2L, tags = c("gene_name","gene_id", "transcript_id","gene_type")))
system.time(gtf_ensembl <- readGFF("~/Documents/HCA/reference/refs/Homo_sapiens.GRCh38.90.gtf", version=2L,tags = c("gene_name","gene_id", "transcript_id","gene_biotype")))
system.time(gtf_refseq <- readGFF("~/Documents/HCA/reference/refs/ncbi-genomes-2017-10-05/GCF_000001405.37_GRCh38.p11_genomic.gff.gz", tags = c("ID", "Name","gbkey","gene","gene_biotype","Parent")))

g4<-unique(gtf_gencode_basic$gene_name)
g3<-unique(gtf_gencode_comp$gene_name)
g2<-unique(gtf_ensembl$gene_name)
g1<-unique(na.omit(gtf_refseq$gene))

## ------------------------------------------------------------------------
ptype<-c('processed_pseudogene','pseudogene','transcribed_unitary_pseudogene','transcribed_unprocessed_pseudogene','unprocessed_pseudogene')
g3.pseudo<-subset(gtf_gencode_comp,gtf_gencode_comp$gene_type %in% ptype & gtf_gencode_comp$type == "gene")
g3.tab<-as.matrix((table(g3.pseudo$gene_type)))
g3.tab

## ------------------------------------------------------------------------
##two mapping categories, uniquely and multiple mapped reads
gene.counts<-list('unq'=list(),'mult'=list())
for(aln in c('unq','mult')){
  ## two read lengths, 25bp and 100bp
  pcounts<-list('25'=c(),'100'=c())
  sralist<-c()
  for(nb in c('25','100')){
    ##load files 
    files<-list.files(path='~/Documents/HCA/reference/counts/',pattern=paste(nb,"_GRCh38_GencodeV27.gene.",aln,".counts.txt",sep=''))
    for(fn in files){
      ## parse out sra ID
      sra<-unlist(strsplit(fn,split='_'))[1]
      sralist<-c(sralist,sra)
      df<-read.delim(paste('~/Documents/HCA/reference/counts/',fn,sep=''),sep='\t',header=T,skip=1)
      ## only need first and seventh column which are gene id, such ensemblID and the actualy counts column.
      x<-df[,c(1,7)]
      colnames(x)<-c('gene_id','counts')
      ## only need the gene annotation, not transcripts
      y<-subset(gtf_gencode_comp,type == 'gene')
      ## left join two table by gene_id
      z<-join(x,y,by='gene_id')
      ## aggregate sum by genetype, such as coding, lncRNA
      z.agg<-aggregate(z$counts,by=list(z$gene_type), FUN=sum,na.rm=T)
      colnames(z.agg)<-c('Group',sra)
      if(length(pcounts[[nb]])==0){
        pcounts[[nb]]<-z.agg
      }else{
        pcounts[[nb]]<-merge(pcounts[[nb]],z.agg,by='Group')
      }
    }
  }
  gene.counts[[aln]]<-pcounts
}

## ------------------------------------------------------------------------
gene.counts.basic<-list('unq'=list(),'mult'=list())
for(aln in c('unq','mult')){
pcounts.basic<-list('25'=c(),'100'=c())
sralist<-c()
for(nb in c('25','100')){
  files<-list.files(path='~/Documents/HCA/reference/counts/',pattern=paste(nb,"_GRCh38_GencodeV27_basic.gene.",aln,".counts.txt",sep=''))
  for(fn in files){
    sra<-unlist(strsplit(fn,split='_'))[1]
    sralist<-c(sralist,sra)
    df<-read.delim(paste('~/Documents/HCA/reference/counts/',fn,sep=''),sep='\t',header=T,skip=1)
    x<-df[,c(1,7)]
    colnames(x)<-c('gene_id','counts')
    y<-subset(gtf_gencode_basic,type == 'gene')
    z<-join(x,y,by='gene_id')
    z.agg<-aggregate(z$counts,by=list(z$gene_type), FUN=sum,na.rm=T)
    colnames(z.agg)<-c('Group',sra)
    if(length(pcounts.basic[[nb]])==0){
      pcounts.basic[[nb]]<-z.agg
    }else{
      pcounts.basic[[nb]]<-merge(pcounts.basic[[nb]],z.agg,by='Group')
    }
}
}
gene.counts.basic[[aln]]<-pcounts.basic
}

## ------------------------------------------------------------------------
atype<-c('protein_coding','antisense_RNA','processed_transcript','lincRNA','Mt_rRNA')
rlen<-c('100')
alns<-c('unq','mult')
for (rr in rlen){
  output<-c()
  output.per<-c()
  for(aln in alns ){
    pcounts<-gene.counts[[aln]]
    pcounts.basic<-gene.counts.basic[[aln]]
    ## load table
    s1<-pcounts[[rr]]
    s2<-pcounts.basic[[rr]]
  ## total gene counts from all annotation
    s1.tot<-apply(s1[,-1],2,sum)
    s2.tot<-apply(s2[,-1],2,sum)
  ## pseudo gene counts
    p1<-subset(s1,s1$Group %in% ptype)
    p2<-subset(s2,s2$Group %in% ptype)
    p1.sum<-data.frame('Group'='pseudo',t(apply(p1[,-1],2,sum)))
    p2.sum<-data.frame('Group'='pseudo',t(apply(p2[,-1],2,sum)))
  ## total counts in each gene category
    c1<-subset(s1,s1$Group %in% atype)
    c2<-subset(s2,s2$Group %in% atype)
  ## others, not in listed annotation types or pseudogene
    o1<-subset(s1,!(s1$Group %in% c(ptype,atype)))
    o2<-subset(s2,!(s2$Group %in% c(ptype,atype)))
    ## sum of category
    o1.sum<-data.frame('Group'='others',t(apply(o1[,-1],2,sum)))
    o2.sum<-data.frame('Group'='others',t(apply(o2[,-1],2,sum)))
  ## total gene counts
    tot1<-rbind(p1.sum,c1,o1.sum)
    tot2<-rbind(p2.sum,c2,o2.sum)
    tot1.ct<-data.frame('Mapped'=rep(aln,nrow(tot1)),'Annot'=rep('comp',nrow(tot1)),tot1)
    tot2.ct<-data.frame('Mapped'=rep(aln,nrow(tot2)),'Annot'=rep('comp',nrow(tot2)),tot2)
    output<-rbind(output,tot1.ct,tot2.ct)
 ## percentage 
    tot1.per<-data.frame('Mapped'=rep(aln,nrow(tot1)),'Annot'=rep('comp',nrow(tot1)),'Group'=tot1[,1],t(apply(tot1[,-1],1,function(x){x/s1.tot})))
    tot2.per<-data.frame('Mapped'=rep(aln,nrow(tot2)),'Annot'=rep('basic',nrow(tot2)),'Group'=tot2[,1],t(apply(tot2[,-1],1,function(x){x/s2.tot})))
    output.per<-rbind(output.per,tot1.per,tot2.per)
  }
  ## visualize
  ##for(g in unique(output.per$Group)){
  ##  x<-subset(output.per,Group == g)
  ##  p<-ggplot(data=melt(x),aes(x=variable,y=value,color=Annot,shape=Mapped))+geom_point(size=4)
  ##  p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
  ##  p<-p+xlab(' ')+ylab(paste('ratio of ',g,' counts in total gene counts',sep=''))+ggtitle(label = paste('Counts of ',g,' type',sep=''))
  ##  p<-p+theme(plot.title = element_text(hjust = 0.5))
  ##  ggsave(p,file=paste('/Users/jishuxu/Documents/HCA/reference/plots/',g,'_rlen_',rr,'_gencode_percent_gene_counts_annotation.png',sep=''),type='cairo-png')
  ##}
  ##write.table(output.per,file=paste('~/Documents/HCA/reference/rlen_',rr,'_gencode_percent_gene_counts_annotation.csv',sep=''),sep=',',col.names=T,row.names=F,quote=F)
  ##write.table(output,file=paste('~/Documents/HCA/reference/rlen_',rr,'_gencode_gene_counts_annotation.csv',sep=''),sep=',',col.names=T,row.names=F,quote=F)
}

## ------------------------------------------------------------------------
gencode.pseudo<-subset(output.per,Group == 'pseudo')
p<-ggplot(data=melt(gencode.pseudo),aes(x=variable,y=value,color=Annot,shape=Mapped))+geom_point(size=4)
p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
p<-p+xlab(' ')+ylab('sum(pseudogene counts)/sum(all gene counts)')+ggtitle(label = paste('% of pseudogene',sep=''))
p<-p+theme(plot.title = element_text(hjust = 0.5))
p

## ------------------------------------------------------------------------
x<-subset(output.per,Group == 'protein_coding')
p<-ggplot(data=melt(x),aes(x=variable,y=value,color=Annot,shape=Mapped))+geom_point(size=4)
p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
p<-p+xlab(' ')+ylab("sum(protein coding counts)/sum(all gene counts)")+ggtitle(label = paste('% of protein coding',sep=''))
p<-p+theme(plot.title = element_text(hjust = 0.5))
p
gencodes.coding<-x

## ------------------------------------------------------------------------
gene.counts.refseq<-list('unq'=list(),'mult'=list())
for(aln in c('unq','mult')){
  pcounts<-list('25'=c(),'100'=c())
  sralist<-c()
  for(nb in c('25','100')){
    files<-list.files(path='~/Documents/HCA/reference/counts/',pattern=paste(nb,"_GRCh38_RefSeq.gene.",aln,".counts.txt",sep=''))
    for(fn in files){
      sra<-unlist(strsplit(fn,split='_'))[1]
      sralist<-c(sralist,sra)
      df<-read.delim(paste('~/Documents/HCA/reference/counts/',fn,sep=''),sep='\t',header=T,skip=1)
      x<-df[,c(1,7)]
      colnames(x)<-c('ID','counts')
      y<-as.data.frame(subset(gtf_refseq,type == 'gene'))
      z<-join(x,y,by='ID')
      z.agg<-aggregate(z$counts,by=list(z$gene_biotype), FUN=sum,na.rm=T)
      colnames(z.agg)<-c('Group',sra)
      if(length(pcounts[[nb]])==0){
        pcounts[[nb]]<-z.agg
      }else{
        pcounts[[nb]]<-merge(pcounts[[nb]],z.agg,by='Group')
      }
    }
  }
  gene.counts.refseq[[aln]]<-pcounts
}

## ------------------------------------------------------------------------
atype<-c('protein_coding','processed_transcript','lncRNA','Mt_rRNA','rRNA','pseudogene','transcribed_pseudogene')
rlen<-c('100')
alns<-c('unq','mult')
for (rr in rlen){
  output<-c()
  output.per<-c()
  for(aln in alns ){
    pcounts<-gene.counts.refseq[[aln]]
    s1<-pcounts[[rr]]
  ## total gene counts
    s1.tot<-apply(s1[,-1],2,sum)
  ## tota protein coding gene counts
    c1<-subset(s1,s1$Group %in% atype)
  ## others
    o1<-subset(s1,!(s1$Group %in% atype))
    o1.sum<-data.frame('Group'='others',t(apply(o1[,-1],2,sum)))
  ## total gene counts
    tot1<-rbind(c1,o1.sum)
    tot1.ct<-data.frame('Mapped'=rep(aln,nrow(tot1)),'Annot'=rep('RefSeq',nrow(tot1)),tot1)
    output<-rbind(output,tot1.ct)
 ## percentage
    tot1.per<-data.frame('Mapped'=rep(aln,nrow(tot1)),'Annot'=rep('RefSeq',nrow(tot1)),'Group'=tot1[,1],t(apply(tot1[,-1],1,function(x){x/s1.tot})))
    output.per<-rbind(output.per,tot1.per)
  }
  ## visualize
  
  ##for(g in unique(output.per$Group)){
  ##  x<-subset(output.per,Group == g)
  ##  p<-ggplot(data=melt(x),aes(x=variable,y=value,color=Annot,shape=Mapped))+geom_point(size=4)
  ##  p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
   ## p<-p+xlab(' ')+ylab(paste('ratio of ',g,' counts in total gene counts',sep=''))+ggtitle(label = paste('Counts of ',g,' type',sep=''))
  ##  p<-p+theme(plot.title = element_text(hjust = 0.5))
  ##  ggsave(p,file=paste('/Users/jishuxu/Documents/HCA/reference/plots/',g,'_rlen_',rr,'_refseq_percent_gene_counts_annotation.png',sep=''),type='cairo-png')
  ##}
  write.table(output.per,file=paste('~/Documents/HCA/reference/rlen_',rr,'_refseq_percent_gene_counts_annotation.csv',sep=''),sep=',',col.names=T,row.names=F,quote=F)
  write.table(output,file=paste('~/Documents/HCA/reference/rlen_',rr,'_refseq_gene_counts_annotation.csv',sep=''),sep=',',col.names=T,row.names=F,quote=F)
}

## ------------------------------------------------------------------------
x<-subset(output.per,Group %in% c('pseudogene','transcribed_pseudogene'))
refseq.pseudo<-aggregate(x[,-c(1,2,3)],by=list(x$Mapped),FUN=sum)
y<-melt(aggregate(x[,-c(1,2,3)],by=list(x$Mapped),FUN=sum))
colnames(y)<-c('Group','Sample','value')
p<-ggplot(data=y,aes(x=Sample,y=value,color=Group,shape=Group))+geom_point(size=4)
p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
p<-p+xlab(' ')+ylab('ratio of pseudo gene counts in total gene counts')+ggtitle(label = paste('Counts of pseudo gene  type',sep=''))
p<-p+theme(plot.title = element_text(hjust = 0.5))
p

## ------------------------------------------------------------------------
refseq.coding<-subset(output.per,Group %in% c('protein_coding') )
y<-melt(aggregate(refseq.coding[,-c(1,2,3)],by=list(refseq.coding$Mapped),FUN=sum))
colnames(y)<-c('Group','Sample','value')
p<-ggplot(data=y,aes(x=Sample,y=value,color=Group,shape=Group))+geom_point(size=4)
p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
p<-p+xlab(' ')+ylab('ratio of coding gene counts in total gene counts')+ggtitle(label = paste('% of coding gene  type',sep=''))
p<-p+theme(plot.title = element_text(hjust = 0.5))
p


## ------------------------------------------------------------------------
xx<-data.frame('Mapped'=refseq.pseudo[,1],'Annot'=rep('RefSeq',2),'Group'=rep('pseudo',2),refseq.pseudo[,2:11])
yy<-rbind(xx,gencode.pseudo)
p<-ggplot(data=melt(yy),aes(x=variable,y=value,color=Annot,shape=Mapped))+geom_point(size=4)
p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
p<-p+xlab(' ')+ylab('sum(pseduogene counts)/sum(all gene counts)')+ggtitle(label = paste('% of pseudogene',sep=''))
p<-p+theme(plot.title = element_text(hjust = 0.5))
p

## ------------------------------------------------------------------------
xx<-refseq.coding
yy<-rbind(gencodes.coding,xx)
p<-ggplot(data=melt(yy),aes(x=variable,y=value,color=Annot,shape=Mapped))+geom_point(size=4)
p<-p+theme(axis.text.y = element_text(color='black',size=10,face='bold'),axis.text.x = element_text(color='black',size=10,face='bold',hjust=1,angle = 45))
p<-p+xlab(' ')+ylab('sum(protein counts)/sum(all gene counts)')+ggtitle(label = paste('% of protein coding',sep=''))
p<-p+theme(plot.title = element_text(hjust = 0.5))
p


