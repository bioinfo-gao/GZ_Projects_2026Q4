suppressMessages({library(scDblFinder);library(SingleCellExperiment);library(Matrix)})
cat("scDblFinder",as.character(packageVersion("scDblFinder")),"xgboost",as.character(packageVersion("xgboost")),"\n")
sim<-function(k,nc=3000,ng=2000,fd=.10,seed=1){set.seed(seed)
 base<-rlnorm(ng,-.5,1.2); mus<-sapply(1:k,function(i){f<-rep(1,ng); if(k>1){j<-sample(ng,150); f[j]<-rlnorm(150,2,.5)}; base*f})
 lab<-sample(k,nc,TRUE); lib<-rlnorm(nc,0,.3)
 M<-sapply(seq_len(nc),function(i) rnbinom(ng,mu=mus[,lab[i]]*lib[i],size=2))
 nd<-round(nc*fd); a<-sample(nc,nd); b<-sample(setdiff(seq_len(nc),a),nd)
 if(k>1){ # 只保留跨类型(heterotypic)配对
   keep<-lab[a]!=lab[b]; a<-a[keep]; b<-b[keep]; nd<-length(a)}
 D<-M[,a]+M[,b]; X<-cbind(M[,-unique(c(a,b))],D); y<-c(rep(0,ncol(X)-nd),rep(1,nd))
 dimnames(X)<-list(paste0("g",1:ng),paste0("c",1:ncol(X))); list(X=as(X,"dgCMatrix"),y=y)}
auprc<-function(s,y){o<-order(-s);tp<-cumsum(y[o]);p<-tp/seq_along(o);sum(p[y[o]==1])/sum(y)}
for(k in c(3,1)) for(seed in 1:2){d<-sim(k,seed=seed)
 s<-scDblFinder(SingleCellExperiment(list(counts=d$X)),BPPARAM=BiocParallel::SerialParam(),verbose=FALSE)
 call<-s$scDblFinder.class=="doublet"; 
 cat(sprintf("k=%d seed=%d n=%d true_dbl=%d called=%d recall=%.3f precision=%.3f AUPRC=%.3f base=%.3f auc_libsize_only=%.3f\n",k,seed,ncol(d$X),sum(d$y),sum(call),mean(call[d$y==1]),ifelse(sum(call)>0,mean(d$y[call]==1),NA),auprc(s$scDblFinder.score,d$y),mean(d$y),auprc(Matrix::colSums(d$X),d$y)))}
