#' transcriptionalVarianceExplained
#'
#' Computes for each gene the variance explained by each variable.
#'
#' @param unit Either "tpm" (default), "fpkm", or "counts".
#' @param logTransform Logical; whether to log-transform the data (default T).
#' @param variables The explanatory variables to include. Should be any combination of "sex", "individual", "clone", and "replicate".
#' @param doPlot Logical; whether to plot a summary of the results (default TRUE; requires packages `beanplot` and `LSD` to be installed)
#' @param ncores The number of threads (by defaults uses all but 1 cores)
#'
#' @return a dataframe containing the square residuals, with genes as rows and explanatory variables as columns (plus a first column containing the gene's mean expression).
#'
#' @export
transcriptionalVarianceExplained <- function(unit="tpm", logTransform=T, variables=c("sex","individual","clone","replicate"), doPlot=TRUE, ncores=NULL){
    variables <- match.arg(variables, c("sex","individual","clone","replicate"), several.ok=T)
    unit <- match.arg(unit, c("tpm","fpkm","counts"))
    data("hipsci_annotation")
    annotation$replicate <- paste(annotation$clone,1:nrow(annotation),sep=".")
    egn <- TMM(getGeneExpr(unit))
    if(logTransform) egn <- log(egn+1)
    d <- data.frame(row.names=row.names(egn), meanExp=rowMeans(egn))
    names(d)[1] <- ifelse(logTransform, paste("mean","log",unit,sep="."), paste("mean",unit,sep="."))
    m <- annotation[,variables]
    for(i in 1:ncol(m)) m[,i] <- as.factor(m[,i])
    
    message("Analyzing variance for each gene...")
    cl <- initPar(ncores)
    if(length(cl)>1){
        clusterExport(cl, c("getAOVresiduals"))
        d <- cbind(d,t(parApply(cl, egn, 1, m=m, FUN=getAOVresiduals)))
        stopCluster(cl)
    }else{
        d <- cbind(d,t(apply(egn,1,m=m,FUN=getAOVresiduals)))
    }
    names(d)[1+1:length(variables)] <- variables
    if(doPlot){
        if(!.checkPkg("beanplot") | !.checkPkg("LSD")){
            warning("Skipping plots: packages 'beanplot' and 'LSD' must be installed for plotting the results...")
        }else{
            layout(matrix(1:2,nrow=2))
            ve <- colSums(d[,-1])
            ve <- ve/sum(ve)
            LSD::LSD.pie(ve,addPercent=T,labels=names(ve),main="")
            tmp <- as.data.frame(d[,-1]/rowSums(d[,-1]))
            beanplot::beanplot(tmp,ll=0.06,col=c("black","black","grey","red"),what=c(0,1,1,0),xlab="Proportion of the variance explained",main="Gene-wise analysis of variance",log="",horizontal=T)
        }
    }
    return(d)
}



#' cellphenoVarianceExplained
#'
#' Computes for each cellular feature the variance explained by each variable.
#'
#' @param logTransform Logical; whether to log-transform the data (default TRUE).
#' @param doPlot Logical; whether to plot a summary of the results (default TRUE; requires package `pheatmap` to be installed)
#'
#' @return a dataframe containing the square residuals, with cellular features as rows and explanatory variables as columns
#'
#' @export
cellphenoVarianceExplained <- function(logTransform=TRUE, doPlot=TRUE){
    data("cellpheno")
    anno <- cellpheno[,1:2]
    anno$replicate <- paste(cellpheno$clone,1:nrow(cellpheno),sep="_")
    anno <- anno[,c(2,1,3)]
    pheno <- cellpheno[,-1:-2]
    if(logTransform) pheno <- log(pheno)
    d <- t(apply(t(pheno),1,m=anno,FUN=getAOVresiduals))
    colnames(d) <- colnames(anno)
    if(doPlot){
        if(!.checkPkg("pheatmap")){
            warning("Skipping plots: package 'pheatmap' must be installed for plotting the results...")
        }else{
	    pheatmap(100*d/rowSums(d),col=colorRampPalette(c("white","red","darkred"))(100),cluster_rows=F,cluster_cols=F,display_numbers=T,number_format="%.0f%%",number_color="black",border_color=NA,main="Variance explained by feature and variable")
        }
    }
    return(d)
}




#' getAOVresiduals
#'
#' Returns the residuals of x explained by each variable (column of m).
#'
#' @param x a numeric vector (expression values to be explained)
#' @param m a data.frame with rows corresponding to each value of x, and columns indicating the values of explanatory variables.
#'
#' @return a numeric vector of squared residuals.
#'
#' @export
getAOVresiduals <- function(x,m){ 
    a <- try(aov(as.numeric(x)~.,data=m),silent=T)
    if("try-error" %in% class(a)) return(rep(NA,ncol(m)))
    summary(a)[[1]][,2]
}