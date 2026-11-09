suppressMessages({library(DESeq2);library(clusterProfiler);library(org.Mm.eg.db);library(dplyr)})
dds<-readRDS("rds_cache/dds.rds"); rl<-readRDS("rds_cache/res_list.rds"); C<-readRDS("rds_cache/contrasts.rds")
run<-function(k){
  r<-results(dds,contrast=c("Group",C[[k]]),independentFiltering=FALSE)
  d<-rl[[k]]; r<-r[d$gene_id,]
  map<-bitr(sub("\\..*","",d$gene_name),"SYMBOL","ENTREZID",org.Mm.eg.db)
  e<-map$ENTREZID[match(d$gene_name,map$SYMBOL)]
  mk<-function(v){ok<-!is.na(e)&is.finite(v); x<-tapply(v[ok],e[ok],function(z) z[which.max(abs(z))]); x<-setNames(as.numeric(x),names(x)); sort(x,decreasing=TRUE)}
  s1<-sign(d$log2FoldChange)*-log10(d$pvalue); s1[!is.finite(s1)]<-sign(d$log2FoldChange[!is.finite(s1)])*330
  out<-c()
  for(nm in c("signlogp","wald")){
    gl<-mk(if(nm=="signlogp") s1 else r$stat); set.seed(42)
    x<-suppressWarnings(gseGO(gl,OrgDb=org.Mm.eg.db,ont="BP",keyType="ENTREZID",pvalueCutoff=1,pAdjustMethod="BH",seed=TRUE,verbose=FALSE))
    df<-as.data.frame(x); out[nm]<-sprintf("%s: sig=%d minpadj=%.2g top=%s",nm,sum(df$p.adjust<0.05),min(df$p.adjust),paste(head(df$Description[order(df$p.adjust)],3),collapse=" | "))
  }
  cat(k,"\n ",out[1],"\n ",out[2],"\n")
}
for(k in c("Crtc2-PGE2_vs_control-PGE2","control-PGE2_vs_control-un")) run(k)
