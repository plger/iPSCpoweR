#' edgeRwrapper
#'
#' Performs one differential expression analysis using edgeR.
#'
#' @param e the expression matrix
#' @param mm the design model (as created by the `model.matrix()' function), for paired analyses only. Default NULL.
#' @param DEsamples the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.
#' @param paired Logical; whether the analysis should be paired or not.
#' @param nested Ignored; for consistency with voomWrapper.
#'
#' @return A data.frame.
#'
#' @export
edgeRwrapper <- function(e, mm=NULL, DEsamples=NULL, paired=FALSE, nested=NULL){
  if(paired & !is.null(mm)){
    row.names(mm) <- colnames(e)
    dds <- DGEList(e)
    dds <- calcNormFactors(dds)
    dds <- estimateDisp(dds, mm)
    dds <- glmFit(dds, mm)
    res <- as.data.frame(topTags(glmLRT(dds,"group2"),nrow(e)))
  }else{
    if(is.null(DEsamples))  stop("DEsamples is null!")
    groups <- rep("B",ncol(e))
    groups[DEsamples] <- "A"
    dds <- DGEList(e,group=groups)
    dds <- calcNormFactors(dds)
    capture.output({dds <- estimateDisp(dds)})
    res <- as.data.frame(topTags(exactTest(dds),nrow(e)))
  }
  res[row.names(e),]
}



#' voomWrapper
#'
#' Performs one differential expression analysis using voom
#'
#' @param e the expression matrix
#' @param mm the design model (as created by the `model.matrix()' function), for paired analyses only. Default NULL.
#' @param nested NULL or character vector of size=ncol(e), indicating the blocking variable (default NULL).
#' @param DEsamples the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.
#' @param paired Ignored; for consistency with edgeRwrapper.
#'
#' @return A data.frame.
#'
#' @export
voomWrapper <- function(e, mm=NULL, nested=NULL, DEsamples=NULL, paired=NULL){
  if(is.null(mm)){
    d <- data.frame(row.names=colnames(e), group=rep("B",ncol(e)), stringsAsFactors=F)
    d$group[DEsamples] <- "A"
    mm <- model.matrix(~group,data=d)
  }
  y <- calcNormFactors(DGEList(e))
  if(!is.null(nested)){
    v <- voom(y,mm)
    corfit <- duplicateCorrelation(v, mm, block=nested)
    v <- voom(y, mm, block=nested, correlation=corfit$consensus.correlation)
    fit <- lmFit(v, mm, block=nested, correlation=corfit$consensus.correlation)
    fit2 <- eBayes(fit)
    res <- topTable(fit2, coef="groupB", number=nrow(e))
  }else{
    v <- voom(y,mm)
    fit <- lmFit(v, mm)
    fit2 <- eBayes(fit)
    res <- topTable(fit2, coef="groupB", number=nrow(e))
  }
  res <- res[row.names(e),]
  return(data.frame(row.names=row.names(res), logFC=res$logFC, PValue=res$P.Value, FDR=res$adj.P.Val))
}


#' voomWrapperSumReps
#'
#' Performs one differential expression analysis using voom
#'
#' @param e the expression matrix
#' @param mm the design model (as created by the `model.matrix()' function), for paired analyses only. Default NULL.
#' @param nested NULL or character vector of size=ncol(e), indicating the blocking variable (default NULL).
#' @param DEsamples the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.
#' @param paired Ignored; for consistency with edgeRwrapper.
#'
#' @return A data.frame.
#'
#' @export
voomWrapperSumReps <- function(e, mm=NULL, nested=NULL, DEsamples=NULL, paired=NULL){
  if(is.null(mm)){
    d <- data.frame(row.names=colnames(e), group=rep("B",ncol(e)), stringsAsFactors=F)
    d$group[DEsamples] <- "A"
    mm <- model.matrix(~group,data=d)
  }
  if(!is.null(nested)){
    e2 <- data.frame(row.names=row.names(e))
    for(s in unique(nested)) e2[[s]] <- rowSums(e[,which(nested==s)])
    e <- e2
    group <- d$group[which(!duplicated(nested))]
    mm <- model.matrix(~group)
  }
  y <- calcNormFactors(DGEList(e))
  v <- voom(y,mm)
  fit <- lmFit(v, mm)
  fit2 <- eBayes(fit)
  res <- topTable(fit2, coef="groupB", number=nrow(e))
  res <- res[row.names(e),]
  return(data.frame(row.names=row.names(res), logFC=res$logFC, PValue=res$P.Value, FDR=res$adj.P.Val))
}


#' dreamWrapper
#'
#' Performs one differential expression analysis using the dream method from the
#' `variancePatition` package (required).
#'
#' @param e the expression matrix
#' @param mm Ignored; for consistency with edgeRwrapper.
#' @param nested A character vector of size=ncol(e), indicating the individuals.
#' @param DEsamples the index of the samples in which the differential expression was added.
#' @param paired Ignored; for consistency with edgeRwrapper.
#'
#' @return A data.frame.
#'
#' @export
dreamWrapper <- function(e, mm=NULL, nested=NULL, DEsamples=NULL, paired=NULL){
  library(variancePartition)
  stopifnot(!is.null(nested))
  d <- data.frame(row.names=colnames(e), group=rep("B",ncol(e)), stringsAsFactors=F)
  d$group[DEsamples] <- "A"
  d$individual <- nested
  form <- ~group+(1|individual)
  dds <- DGEList(e)
  dds <- calcNormFactors(dds)
  vobjDream = voomWithDreamWeights( dds, form, d )
  fitmm = suppressWarnings(dream( vobjDream, form, d ))
  res <- as.data.frame(topTable( fitmm, coef='groupB', number=Inf ))[row.names(e),]
  return(data.frame(row.names=row.names(res), logFC=res$logFC, PValue=res$P.Value, FDR=res$adj.P.Val))
}



#' voomLmerWrapper
#'
#' Performs one differential expression analysis using mixed models (via `lme4::lmer`) on voom-transformed counts
#'
#' @param e the expression matrix
#' @param nested NULL or character vector of size=ncol(e), indicating the blocking variable (default NULL).
#' @param DEsamples the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.
#' @param mm Ignored; for consistency with alternative functions.
#' @param paired Ignored; for consistency with alternative functions.
#'
#' @return A data.frame.
#'
#' @export
voomLmerWrapper <- function(e, nested=NULL, DEsamples=NULL, mm=NULL, paired=NULL){
  library("lme4")
  group <- rep("B",ncol(e))
  group[DEsamples] <- "A"
  mm <- model.matrix(~group)
  v <- voom(calcNormFactors(DGEList(e)), mm)
  res <- topTable(eBayes(lmFit(v,mm)), coef="groupB", number=nrow(e))
  res <- res[row.names(e),]
  res[,"P.Value"] <- apply(v,1,nested=nested,group=group,FUN=function(x,group,nested){
    x <- as.numeric(x)
    if(is.null(nested)){
      mod <- lm(x~1+group)
      summary(mod)[[4]]["groupB","Pr(>|t|)"]
    }else{
      mod <- try(lmer(x~1+(1|nested)+group),silent=T)
      if(!is(mod,"try-error")){
        a <- drop1(mod, test="Chisq")["group","Pr(Chi)"]
        if(a>=0) return(a)
      }
    }
    return(1)
  })
  res$FDR <- p.adjust(res$P.Value, method="fdr")
  return(data.frame(row.names=row.names(res), logFC=res$logFC, PValue=res$P.Value, FDR=res$FDR))
}


#' vstLmerWrapper
#'
#' Performs one differential expression analysis using mixed models (via `lme4::lmer`) on `DESeq2` vst-transformed counts
#'
#' @param e the expression matrix
#' @param nested NULL or character vector of size=ncol(e), indicating the blocking variable (default NULL).
#' @param DEsamples the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.
#' @param mm Ignored; for consistency with alternative functions.
#' @param paired Ignored; for consistency with alternative functions.
#'
#' @return A data.frame.
#'
#' @export
vstLmerWrapper <- function(e, nested=NULL, DEsamples=NULL, mm=NULL, paired=NULL){
  library(DESeq2)
  library("lme4")
  group <- rep("B",ncol(e))
  group[DEsamples] <- "A"
  for(i in 1:ncol(e)) e[,i] <- as.integer(floor(e[,i]))
  cds <- DESeqDataSetFromMatrix(countData=e, colData=data.frame(row.names=colnames(e), group=as.factor(group)), design=~group)
  cds <- estimateSizeFactors(cds)
  e <- assay(vst(cds, blind=T))
  res <- data.frame(row.names=row.names(e), logFC=rep(NA, nrow(e)))
  res$PValue <- apply(e,1,group=group,nested=nested,FUN=function(x,group,nested){
    x <- as.numeric(x)
    if(is.null(nested)){
      mod <- lm(x~1+group)
      summary(mod)[[4]]["groupB","Pr(>|t|)"]
    }else{
      mod <- try(lmer(x~1+(1|nested)+group),silent=T)
      if(!is(mod,"try-error")){
        a <- drop1(mod, test="Chisq")["group","Pr(Chi)"]
        if(a>=0) return(a)
      }
    }
    return(1)
  })
  res$FDR <- p.adjust(res$PValue, method="fdr")
  return(res)
}


#' glmmWrapper
#'
#' Performs one differential expression analysis using mixed models (via `MASS::glmmPQL`) and quasipoisson distribution.
#'
#' @param e the expression matrix
#' @param nested character vector of size=ncol(e), indicating the blocking variable
#' @param DEsamples the index of the samples in which the differential expression was added. Mandatory if mm=NULL, ignored otherwise.
#' @param mm Ignored; for consistency with alternative functions.
#' @param paired Ignored; for consistency with alternative functions.
#'
#' @return A data.frame.
#'
#' @export
glmmWrapper <- function(e, nested, DEsamples=NULL, mm=NULL, paired=NULL){
  library("MASS")
  group <- rep("B",ncol(e))
  group[DEsamples] <- "A"
  mm <- model.matrix(~group)
  nf <- calcNormFactors(e)
  v <- voom(calcNormFactors(DGEList(e)), mm)
  res <- topTable(eBayes(lmFit(v,mm)), coef="groupB", number=nrow(e))
  res <- res[row.names(e),]
  res[,"P.Value"] <- apply(e,1,nf=nf, group=group, nested=nested,FUN=function(x,nf,group,nested){
    x <- log(as.numeric(x)+1)
    mod <- try(glmmPQL(x~1+nf+group, ~1|nested, family="quasipoisson"),silent=T)
    if(!is(mod,"try-error")){
      return(summary(mod)$tTable[3,5])
    }
    return(1)
  })
  res$FDR <- p.adjust(res$P.Value, method="fdr")
  return(data.frame(row.names=row.names(res), logFC=res$logFC, PValue=res$P.Value, FDR=res$FDR))
}
