suppressMessages({library(Seurat)}); cat("Seurat",as.character(packageVersion("Seurat")),"\n")
set.seed(1); n<-300; g<-60
m<-matrix(rnbinom(n*g,2,.5),g,n,dimnames=list(c(paste0("MT-",1:13),paste0("G",1:47)),paste0("c",1:n)))
m[1:13,1:150]<-m[1:13,1:150]*20L   # 前150细胞（样本A）mt 高
x<-as(m,"dgCMatrix")
truth<-Matrix::colSums(x[1:13,])/Matrix::colSums(x)*100
o<-CreateSeuratObject(x); o$s<-rep(c("A","B"),each=150)
a<-PercentageFeatureSet(o,pattern="^MT-")
o2<-o; o2[["RNA"]]<-split(o2[["RNA"]],f=o2$s); print(Layers(o2))
b<-tryCatch(PercentageFeatureSet(o2,pattern="^MT-"),error=function(e){cat("ERR split:",conditionMessage(e),"\n");NULL})

va<-if(is.data.frame(a)) setNames(a[,1],rownames(a)) else a
cat("class a:",class(a),"\n")
cat("unsplit max|diff|:",max(abs(va[names(truth)]-truth)),"\n")
if(!is.null(b)){vb<-if(is.data.frame(b)) setNames(b[,1],rownames(b)) else b; cat("split max|diff|:",max(abs(vb[names(truth)]-truth),na.rm=TRUE),"NA:",sum(is.na(vb)),"\n")}
cnt<-LayerData(o2,layer="counts.A"); cat("layers ok\n")
suppressMessages(library(scDblFinder)); suppressMessages(library(SingleCellExperiment))
set.seed(2); M<-matrix(rnbinom(2000*1000,mu=1,size=1),1000,2000,dimnames=list(paste0("g",1:1000),paste0("c",1:2000)))
r<-tryCatch({s<-scDblFinder(SingleCellExperiment(list(counts=as(M,"dgCMatrix"))),BPPARAM=BiocParallel::SerialParam(),verbose=FALSE);table(s$scDblFinder.class)},error=function(e)paste("ERR",conditionMessage(e))); print(r)
