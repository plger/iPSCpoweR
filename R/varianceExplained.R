#' transcriptionalVarianceExplained
#'
#' Computes for each gene the variance explained by each variable.
#'
#' @param unit Either "tpm" (default), "fpkm", or "counts".
#' @param logTransform Logical; whether to log-transform the data (default T).
#' @param variables The explanatory variables to include. Default "sex" and "individual"
#' @param doPlot Logical; whether to plot a summary of the results (default TRUE; requires packages `beanplot` and `LSD` to be installed)
#' @param mixed Logical; whether to use mixed models (via `lme4`) (default TRUE).
#' @param formula Specify the formula to use. If specified, this overrides the `variables` argument. If `mixed=F`, the default formula is `~.`; if `mixed=T`, each variable is considered a random variable.
#' @param ncores The number of threads (by defaults uses all but 1 cores)
#'
#' @return a dataframe containing the square residuals, with genes as rows and explanatory variables as columns (plus a first column containing the gene's mean expression).
#'
#' @export
transcriptionalVarianceExplained <- function(res=NULL, unit="tpm", logTransform=T, variables=c("sex","individual"), doPlot=TRUE, mixed=TRUE, formula=NULL, ncores=NULL){
    unit <- match.arg(unit, c("tpm","fpkm","counts"))
    if(is.null(res)){
	data("hipsci_annotation")
	annotation <- annotation[which(!annotation$feeder),]
	egn <- TMM(getGeneExpr(unit))
	egn <- egn[,row.names(annotation)]
    }else{
	annotation <- res$annotation
	egn <- res$dat
	rm(res)
    }
    if(logTransform) egn <- log(egn+1)
    colnames(annotation) <- tolower(colnames(annotation))
    variables <- tolower(variables)
    if(!all(variables %in% colnames(annotation))){
	missing <- variables[which(!(variables %in% colnames(annotation)))]
	stop(paste("The following requested variables could not be found in the annotation:",paste(missing,collapse=", ")))
    }
    d <- data.frame(row.names=row.names(egn), meanExp=rowMeans(egn))
    names(d)[1] <- ifelse(logTransform, paste("mean","log",unit,sep="."), paste("mean",unit,sep="."))
    m <- annotation[,variables]
    for(i in 1:ncol(m)) m[,i] <- as.factor(m[,i])

    if(mixed & is.null(formula)){
	formula <- paste0("~1+",paste0("(1|",variables,")",collapse="+"))
    }
    message("Analyzing variance for each gene...")
    cl <- initPar(ncores)
    if(length(cl)>1){
        clusterExport(cl, c("getLmeResiduals","getAOVresiduals"))
        if(mixed){
		clusterEvalQ(cl, library(lme4))
		tmp <- t(parApply(cl, egn, 1, m=m, variables=variables, form=formula, FUN=getLmeResiduals))
	}else{
		tmp <- t(parApply(cl, egn, 1, m=m, form=formula, FUN=getAOVresiduals))
	}
        stopCluster(cl)
    }else{
        if(mixed){
		library(lme4)
		tmp <- t(apply(egn,1,m=m,variables=variables,form=formula,FUN=getLmeResiduals))
	}else{
		tmp <- t(apply(egn,1,m=m,form=formula,FUN=getAOVresiduals))
	}
    }
    colnames(tmp) <- c(variables,"residuals")
    d <- cbind(d,tmp)
    rm(tmp)
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
#' Returns the residuals of x explained by each variable (column of m) using classical analysis of variance.
#'
#' @param x a numeric vector (expression values to be explained)
#' @param m a data.frame with rows corresponding to each value of x, and columns indicating the values of explanatory variables.
#' @param form Specify the formula to use (default NULL for ~.)
#'
#' @return a numeric vector of squared residuals.
#'
#' @export
getAOVresiduals <- function(x,m,form=NULL){ 
    if(is.null(form)){
	a <- try(aov(as.numeric(x)~.,data=m),silent=T)
    }else{
	a <- try(aov(as.formula(paste("as.numeric(x)",form,sep="")),data=m),silent=T)
    }
    if("try-error" %in% class(a)) return(rep(NA,ncol(m)+1))
    summary(a)[[1]][,2]
}


#' getLmeResiduals
#'
#' Returns the residuals of x explained by each variable (column of m) using mixed models.
#'
#' @param x a numeric vector (expression values to be explained)
#' @param m a data.frame with rows corresponding to each value of x, and columns indicating the values of explanatory variables.
#' @param form Specify the formula to use (required).
#' @param variables Character vector containing the names of the variables, in the order in which they should be returned. If NULL, `colnames(m)` is used.
#'
#' @return a numeric vector of squared residuals.
#'
#' @export
getLmeResiduals <- function(x,m,form,variables=NULL){
    library(lme4)
    if(is.null(variables)) variables <- colnames(m)
    m <- cbind(as.numeric(x),m)
    colnames(m)[1] <- "x"
    form <- as.formula(paste0("x",form))
    mod <- lmer(form,data=m)
    vc <- as.data.frame(VarCorr(mod))
    row.names(vc) <- vc[,1]
    rn <- c(variables,"Residual")
    return(vc[c(variables,"Residual"),"vcov"])
}