#' DEA.permutateIndividuals
#'
#' Performs differential-expression analyses between samples from random individuals.
#'
#' @param nbIndividuals Either a positive integer, or a vector of positive integers of length 2, indicating the number of individuals in each group (default 2).
#' @param maxTests The maximum number of tests to run (default 50).
#' @param nbClone The number of iPSC clones to consider from each individual (default 2)
#' @param res The result of `aggByClone(getGeneExpr())`, optional and used only to speed up multiple calls or use a custom dataset.
#' @param useOnlyMulticlone Logical; whether to use only the samples that have more than one clone (default TRUE). Disabling this greatly increase the number of possible combinations, but will make the search for viable combinations _very_ slow. It should be used only if insufficient viable combinations are available (e.g. for nbIndividuals > 5).
#' @param seed The seed (for randomization)
#' @param addDE Logical; whether to add random differentially-expressed genes (default FALSE).
#' @param addDE.foldchanges A numeric vector of foldchanges to add (will use the given foldchanges and their inverse, applied on the first group). Default c(1.1,1.25,1.5,2,3,5).
#' @param addDE.nbPerFC The number of DEGs for each foldchange (default 5).
#' @param doSave Whether to save the results (default TRUE).
#' @param returnResults Whether to return the results (default FALSE).
#' @param ncores Integer; number of cores to use. Defaults to detected number of cores minus one.
#' @param quiet Logical; if TRUE, suppresses progress report (default FALSE)
#' @param filter Function; the filter for genes to be tested, to be applied to all rows of the expression matrix (default: no filter). An example value for the filter argument would be `sum(x>10)>2', which would only include genes that have more than 10 counts in more than 2 samples.
#' @param nested Logical; whether to run the analysis using duplicateCorrelation (default FALSE). Requires `nbClone=2` and `DEAfunc=voomWrapper`.
#' @param fasterCombinations Logical; whether to find combinations rapidly by selecting a pre-defined sex ratio (default FALSE). Useful for large datasets.
#' @param DEAfunc Function; the function used to run the differential expression analysis (by default, the `edgeRwrapper' function of this package).
#'
#' @return Nothing (results save to files), or a data.frame if returnResults=TRUE.
#'
#' @export
DEA.permutateIndividuals <- function(	nbIndividuals=2, 
					maxTests=50, 
					nbClone=2, 
					res=NULL, 
					useOnlyMulticlone=TRUE, 
					seed=1, 
					addDE=FALSE, 
					addDE.foldchanges=c(1.25,1.5,2,3,5), 
					addDE.nbPerFC=10, 
					doSave=TRUE, 
					returnResults=FALSE, 
					ncores=NULL, quiet=FALSE, 
					filter=NULL, 
					nested=FALSE, 
					fasterCombinations=FALSE, 
					DEAfunc=edgeRwrapper){
    library(limma)
    library(edgeR)
    if(length(nbClone) != 1 | !(nbClone %in% c(1,2)))  stop("nbClone should be either 1 or 2.")
    if(nbClone==2 & !useOnlyMulticlone) stop("useOnlyMulticlone cannot be FALSE when nbClone=2")
    if(length(nbIndividuals)==1) nbIndividuals <- c(nbIndividuals,nbIndividuals)
    if(sum(nbIndividuals)>14 && length(unique(nbIndividuals))>1) stop("Unbalanced groups are only implemented for sum of group sizes <=14")
    
    saveToFilePrefix <- paste(nbIndividuals[1],"indiv.vs.",nbIndividuals[2],"indiv.",nbClone,ifelse(nested,".nested",""),sep="")
    titlePrefix <- paste(nbIndividuals[1]," individuals vs ",nbIndividuals[2]," individuals (",ifelse(nbClone==2,"2 clones","1 clone"),"/individual)",sep="")

    DEsamples <- 1:(nbClone*nbIndividuals[1])
    
    if(is.null(res))    res <- aggByClone(getGeneExpr(quiet=quiet))
    m <- res$annotation
    agc <- res$dat
    rm(res)
    if(!is.null(filter)){
	if(!("function" %in% class(filter)))	warning("The filter argument is not a function - will be ignored.")
	agc <- agc[which(apply(agc,1,FUN=filter)),]
    }
    
    set.seed(seed)
    
    if(!quiet)  message(titlePrefix)

    if(nested){
	if(nbClone < 2) stop("nested analysis can only be enabled if nbClone>1")
	nested <- m$individual
	names(nested) <- row.names(m)
    }else{
	nested <- NULL
    }
    
    # select input differentially-expressed genes
    agc <- agc[order(rowMeans(agc)),]
    DEgenes <- .randomDEgenes(nrow(agc), addDE.foldchanges, addDE.nbPerFC)
    DEfcs <- rep(c(addDE.foldchanges,1/addDE.foldchanges),addDE.nbPerFC)
    # select combinations of samples
    if(!("sex" %in% colnames(m))){
	if("Sex" %in% colnames(m)){
		m$sex <- m$Sex
	}else{
		m$sex <- "Male"
	}
    }
    id <- m$sex[!duplicated(m$individual)]
    names(id) <- m$individual[!duplicated(m$individual)]
    t <- table(m$individual)
    if(!quiet) message("Looking for combinations...")
    if( nbClone > 1 | useOnlyMulticlone) t <- t[which(t>1)]
    if( nbClone > 1 | useOnlyMulticlone){
        c2 <- as.data.frame(t(combn(names(t),sum(nbIndividuals))),stringsAsFactors=F)
        if(nrow(c2) > maxTests*20) c2 <- c2[sample(nrow(c2),20*maxTests),]
        c2$sex.balanced <- apply(c2,1,id=id,nbIndividuals=nbIndividuals, FUN=function(x,nbIndividuals,id){ .getSexRatio(x[1:nbIndividuals[1]],id)==.getSexRatio(x[nbIndividuals[1]+1:nbIndividuals[2]],id) })
        c2 <- c2[which(c2$sex.balanced),]
    }else if(sum(nbIndividuals) <= min(sexes <- table(id))){
  	  c2 <- do.call(rbind,lapply(names(sexes),m=m,maxTests=maxTests,nbIndividuals=nbIndividuals,st=sexes,FUN=function(x,m,maxTests,nbIndividuals,st){
    		nn <- ceiling(1.2*maxTests*st[x]/sum(st))
    		c3 <- as.data.frame(t(combn(unique(m$individual[which(m$sex==x)]),sum(nbIndividuals))),stringsAsFactors=F)
    		if(nrow(c3) > (nn)) c3 <- c3[sample(nrow(c3),nn),]
    		return(c3)
    	}))
    }else{
      nbFemales <- sum(nbIndividuals)*(max(sexes)/sum(sexes))
      if((ceiling(nbFemales) %% 2)==0){
        nbFemales <- ceiling(nbFemales)
      }else{
        nbFemales <- floor(nbFemales)
      }
      fem <- rev(names(sort(sexes)))[1]
      nbMales <- sum(nbIndividuals)-nbFemales
      a <- t(sapply( 1:maxTests, u=unique(m$individual[which(m$sex==fem)]), FUN=function(x,u) sort(sample(u, nbFemales))))
      a <- a[!duplicated(cbind(t(apply(a[,seq_len(nbFemales/2)],1,sort)), t(apply(a[,-seq_len(nbFemales/2)],1,sort)))),]
      b <- t(sapply( seq_len(nrow(a)), u=unique(m$individual[which(m$sex!=fem)]), FUN=function(x,u) sort(sample(u, nbMales))))
      c2 <- cbind(a[,seq_len(nbFemales/2)], b[,seq_len(floor(nbMales/2))], 
                  a[,-seq_len(nbFemales/2)], b[,-seq_len(floor(nbMales/2))] )
      c2 <- as.data.frame(c2)
    }
    if(nrow(c2) > maxTests)	c2 <- c2[sample(nrow(c2),maxTests),]

    message(paste("Found ",nrow(c2),"comparisons."))
    
    if(nrow(c2)==0){
        message(paste("No sex-balanced comparison for",paste(nbIndividuals,collapse="vs"),"individuals at",nbClone," clone per individual."))
        if(useOnlyMulticlone) message("Try with useOnlyMulticlone=FALSE")
        return(NULL)
    }
    row.names(c2) <- 1:nrow(c2)
    if(nbClone==2){
        c2$clones <- apply(c2[,1:sum(nbIndividuals)],1,FUN=function(x){ x <- as.character(x); paste(c(row.names(m)[which(m$individual %in% x[1:nbIndividuals[1]])],row.names(m)[which(m$individual %in% x[nbIndividuals[1]+1:nbIndividuals[2]])]),collapse=",")})
    }else{
        c2$clones <- apply(c2[,1:sum(nbIndividuals)],1,FUN=function(x){
            paste(sapply(as.character(x),FUN=function(y){ row.names(m)[sample(which(m$individual == y),1)] }), collapse=",")
        })    
    }
    
    if(!quiet)  message(paste("Running",nrow(c2),"differential expression analyses"))
    
    d.logFC <- as.data.frame(row.names=row.names(agc),matrix(0,nrow=nrow(agc),ncol=nrow(c2)))
    d.PValue <- as.data.frame(row.names=row.names(agc),matrix(1,nrow=nrow(agc),ncol=nrow(c2)))
    d.FDR <- as.data.frame(row.names=row.names(agc),matrix(1,nrow=nrow(agc),ncol=nrow(c2)))
        
    cl <- initPar(ncores)
    if(any(class(cl) != "logical")){
        if(!quiet){
            f <- tempfile()
            writeBin(0,f)
            progBar <- txtProgressBar(0, nrow(c2), file=stdout(), style=3)
            setTxtProgressBar(progBar, 0)
        }else{
            f <- NULL
            progBar <- NULL
        }
        clusterEvalQ(cl, library(edgeR))
        clusterExport(cl, c("agc"), envir=environment())
        d <- parLapply(cl, as.character(c2$clones), DEgenes=DEgenes, nested=nested, DEfcs=DEfcs, DEAfunc=DEAfunc, DEsamples=DEsamples, progBar=progBar, progFile=f, fun=.doDEA)
        stopCluster(cl)
        if(!quiet){
            close(progBar)
            cat("\nDone! Aggregating results...\n")
        }
        for(i in 1:nrow(c2)){
            d.logFC[,i] <- d[[i]]$logFC
            d.PValue[,i] <- d[[i]]$PValue
            d.FDR[,i] <- d[[i]]$FDR            
        }
        rm(d)
    }else{
        for(i in 1:nrow(c2)){
            if(!quiet)  cat(paste(i,"/",nrow(c2),"\n",sep=""))
            x <- strsplit(as.character(c2$clones[i]),",",fixed=T)[[1]]
            res <- .doDEA(agc[,x], DEAfunc=DEAfunc, DEgenes=DEgenes, DEfcs=DEfcs, DEsamples=DEsamples, nested=nested)
            d.logFC[,i] <- res$logFC
            d.PValue[,i] <- res$PValue
            d.FDR[,i] <- res$FDR
        }
    }

    res <- data.frame(row.names=rownames(d.FDR),meanPval=rowMeans(d.PValue), pval.below.05=apply(d.PValue,1,FUN=function(x){ sum(x < 0.05, na.rm=T)}), pval.below.01=apply(d.PValue,1,FUN=function(x){ sum(x < 0.01, na.rm=T)}), FDR.below.05=apply(d.FDR,1,FUN=function(x){ sum(x < 0.05, na.rm=T)}), FDR.below.01=apply(d.FDR,1,FUN=function(x){ sum(x < 0.01, na.rm=T)}))
    if(addDE){
        res$inputDE.foldchange <- 1
        res$inputDE.foldchange[DEgenes] <- DEfcs
    }
    agcn <- TMM(agc)
    res$meanCount <- rowMeans(agcn)
    res$medianCount <- apply(agcn,1,FUN=median)
    res$CV <- apply(agcn,1,FUN=function(x){ x <- as.numeric(x); sd(x)/mean(x)})
    res$median.logFC.below.p05 <- sapply(1:nrow(res),FUN=function(x){ median(abs(as.numeric(d.logFC[x,which(d.PValue[x,]<0.05)])),na.rm=T)})
    res$mean.logFC.below.p05 <- sapply(1:nrow(res),FUN=function(x){ mean(abs(as.numeric(d.logFC[x,which(d.PValue[x,]<0.05)])),na.rm=T)})
    res <- res[order(res$FDR.below.05,decreasing=T),]

    permres <- list(comparisons=c2, PValues=d.PValue, FDR=d.FDR, logFC=d.logFC, byGene=res, seed=seed, call=match.call())

    nsig <- apply(d.FDR,2,FUN=function(x){ sum(x<0.05, na.rm=T)})
    
    if(doSave){
        save(permres, file=paste(saveToFilePrefix,"RData",sep="."))
        if(!quiet) message(paste("Results saved in ",saveToFilePrefix,".RData",sep=""))
        svg(paste(saveToFilePrefix,"Nsig","svg",sep="."),width=5,height=4)
	hist(nsig,breaks=40,xlab="Differentially expressed genes (FDR<0.05)",ylab="Number of comparisons",col="grey",main=titlePrefix)
	abline(v=mean(nsig),lty="dashed",lwd=2,col="yellow")                                         
	legend("topright",legend=c(paste("Mean:",round(mean(nsig)),"DEGs"),paste("Median:",round(median(nsig)),"DEGs")),bty="n")
        dev.off()
        svg(paste(saveToFilePrefix,"logFC","svg",sep="."),width=5,height=4)
	hist(res$mean.logFC.below.p05,breaks=40,xlab="average log2(foldchange) at FDR<0.05",ylab="Number of comparisons",col="grey",main=titlePrefix)
        dev.off()
	if(addDE){
		svg(paste(saveToFilePrefix,"sensitivity","svg",sep="."),width=5,height=5)
		getSensitivityMatrix(readPermResults(permres))
		dev.off()		
	}
    }
    if(returnResults) return(permres)
}

#' DEA.permutateClones
#'
#' Performs differential-expression analyses between groups constituted by corresponding iPSC clones from random individuals (for each individual, one clone is randomly assigned to one group, and the other to the second group).
#'
#' @param nbIndividuals A positive integer indicating the number of individuals to consider (default 2).
#' @param maxTests The maximum number of tests to run (default 50).
#' @param res The result of `aggByClone(getGeneExpr())`, optional and used only to speed up multiple calls or use a custom dataset.
#' @param seed The seed (for randomization)
#' @param addDE Logical; whether to add random differentially-expressed genes (default FALSE).
#' @param addDE.foldchanges A numeric vector of foldchanges to add (will use the given foldchanges and their inverse). Default c(1.1,1.25,1.5,2,3,5).
#' @param addDE.nbPerFC The number of DEGs for each foldchange (default 5).
#' @param paired Logical; whether to run the analysis in a paired way, i.e. `~individual+group` (default FALSE). 
#' @param doSave Whether to save the results (default TRUE).
#' @param returnResults Whether to return the results (default FALSE).
#' @param ncores Integer; number of cores to use. Defaults to detected number of cores minus one.
#' @param quiet Logical; if TRUE, suppresses progress report (default FALSE)
#' @param filter Function; the filter for genes to be tested, to be applied to all rows of the expression matrix (default: no filter). An example value for the filter argument would be `sum(x>10)>2', which would only include genes that have more than 10 counts in more than 2 samples.
#' @param useExactTest Logical; whether to use edgeR's exact test for unpaired analysis. Works only with DEAfunc=edgeRwrapper.
#' @param DEAfunc Function; the function used to run the differential expression analysis (by default, the `edgeRwrapper' function of this package).
#'
#' @return Nothing (results save to files), or a data.frame if returnResults=TRUE.
#'
#' @export
DEA.permutateClones <- function(nbIndividuals=2, maxTests=50, res=NULL, seed=1, addDE=F, addDE.foldchanges=c(1.25,1.5,2,3,5), addDE.nbPerFC=10, paired=F, doSave=T, returnResults=F, ncores=NULL,quiet=F,filter=NULL,useExactTest=T,DEAfunc=edgeRwrapper){
    library(limma)
    library(edgeR)
    
    if(length(nbIndividuals)>1) stop("nbIndividuals should be a single positive integer.")
    
    saveToFilePrefix <- paste("clones.",nbIndividuals,"indiv",ifelse(paired,".paired",""),sep="")
    titlePrefix <- paste(ifelse(paired,"Paired ",""),"DEA between corresponding clones of ",nbIndividuals," individuals",sep="")
    
    if(is.null(res))	res <- aggByClone(getGeneExpr(quiet=quiet))
    m <- res$annotation
    agc <- res$dat
    rm(res)
    if(!is.null(filter)){
	if(!("function" %in% class(filter)))	warning("The filter argument is not a function - will be ignored.")
	agc <- agc[which(apply(agc,1,FUN=filter)),]
    }

    if(!quiet)  message(titlePrefix)
    
    set.seed(seed)
    
    # select input differentially-expressed genes
    agc <- agc[order(rowMeans(agc)),]
    DEgenes <- .randomDEgenes(nrow(agc), addDE.foldchanges, addDE.nbPerFC)
    DEfcs <- rep(c(addDE.foldchanges,rev(1/addDE.foldchanges)),addDE.nbPerFC)
    
    # select combinations of samples
    t <- table(m$individual)
    c2 <- as.data.frame(t((combn(names(t[which(t>1)]),nbIndividuals))),stringsAsFactors=F)
    set.seed(seed)
    if(nrow(c2) > maxTests){
        c2 <- c2[sample(nrow(c2),maxTests),]
    }
    if(nrow(c2)==0){
        message(paste("No comparison for",nbIndividuals,"individuals."))
        return(NULL)
    }
    
    row.names(c2) <- 1:nrow(c2)
    c2$clones <- apply(c2[,1:nbIndividuals],1,FUN=function(x){
        paste(as.character(t(sapply(as.character(x),FUN=function(y){ sample(row.names(m)[which(m$individual==y)],2) }))),collapse=",")
        })
        
    if(paired){
        indiv <- as.factor(rep(1:nbIndividuals,2))
        group <- as.factor(rep(1:2,each=nbIndividuals))
	mm <- model.matrix(~indiv+group)
	if(!quiet)  message("Using paired analysis with the following design:")
	if(!quiet)  print(mm)
    }else{
        if(useExactTest){
		mm <- NULL
	}else{
		mm <- model.matrix(~group)
	}
    }
    if(!quiet)  message(paste("Running",nrow(c2),"differential expression analyses..."))
    
    d.logFC <- as.data.frame(row.names=row.names(agc),matrix(0,nrow=nrow(agc),ncol=nrow(c2)))
    d.PValue <- as.data.frame(row.names=row.names(agc),matrix(1,nrow=nrow(agc),ncol=nrow(c2)))
    d.FDR <- as.data.frame(row.names=row.names(agc),matrix(1,nrow=nrow(agc),ncol=nrow(c2)))

    cl <- initPar(ncores)
    if(any(class(cl) != "logical")){
        if(!quiet){
            f <- tempfile()
            writeBin(0,f)
            progBar <- txtProgressBar(0, nrow(c2), file=stdout(), style=3)
            setTxtProgressBar(progBar, 0)
        }else{
            f <- NULL
            progBar <- NULL
        }
        clusterEvalQ(cl, library(edgeR))
        clusterExport(cl, c("agc"), envir=environment())
        d <- parLapply(cl, as.character(c2$clones), mm=mm, DEgenes=DEgenes, DEfcs=DEfcs, DEAfunc=DEAfunc, DEsamples=1:nbIndividuals[1], paired=paired, nested=NULL, progBar=progBar, progFile=f, fun=.doDEA)
        stopCluster(cl)
        if(!quiet){
            close(progBar)
            cat("\nDone! Aggregating results...\n")
        }
        for(i in 1:nrow(c2)){
            d.logFC[,i] <- d[[i]]$logFC
            d.PValue[,i] <- d[[i]]$PValue
            d.FDR[,i] <- d[[i]]$FDR            
        }
        rm(d)
    }else{
        for(i in 1:nrow(c2)){
            if(!quiet)  cat(paste(i,"/",nrow(c2),"\n",sep=""))
            x <- strsplit(as.character(c2$clones[i]),",",fixed=T)[[1]]
            res <- .doDEA(agc[,x], mm=mm, DEgenes=DEgenes, DEfcs=DEfcs, DEsamples=1:nbIndividuals[1], paired=paired, nested=NULL, DEAfunc=DEAfunc)
            d.logFC[,i] <- res$logFC
            d.PValue[,i] <- res$PValue
            d.FDR[,i] <- res$FDR
        }
    }

    res <- data.frame(row.names=rownames(d.FDR),meanPval=rowMeans(d.PValue), pval.below.05=apply(d.PValue,1,FUN=function(x){ sum(x < 0.05, na.rm=T)}), pval.below.01=apply(d.PValue,1,FUN=function(x){ sum(x < 0.01, na.rm=T)}), FDR.below.05=apply(d.FDR,1,FUN=function(x){ sum(x < 0.05, na.rm=T)}), FDR.below.01=apply(d.FDR,1,FUN=function(x){ sum(x < 0.01, na.rm=T)}))

    if(addDE){
        res$inputDE.foldchange <- 1
        res$inputDE.foldchange[DEgenes] <- DEfcs
    }

    agcn <- TMM(agc)
    res$meanCount <- rowMeans(agcn)
    res$medianCount <- apply(agcn,1,FUN=median)
    res$CV <- apply(agcn,1,FUN=function(x){ x <- as.numeric(x); sd(x)/mean(x)})
    res$median.logFC.below.p05 <- sapply(1:nrow(res),FUN=function(x){ median(abs(as.numeric(d.logFC[x,which(d.PValue[x,]<0.05)])),na.rm=T)})
    res$mean.logFC.below.p05 <- sapply(1:nrow(res),FUN=function(x){ mean(abs(as.numeric(d.logFC[x,which(d.PValue[x,]<0.05)])),na.rm=T)})
    res <- res[order(res$FDR.below.05,decreasing=T),]

    permres <- list(comparisons=c2, PValues=d.PValue, FDR=d.FDR, logFC=d.logFC, byGene=res, seed=seed, call=match.call())

    nsig <- apply(d.FDR,2,FUN=function(x){ sum(x<0.05, na.rm=T)})
    
    if(doSave){
        save(permres, file=paste(saveToFilePrefix,"RData",sep="."))
        if(!quiet) message(paste("Results saved in ",saveToFilePrefix,".RData",sep=""))
        svg(paste(saveToFilePrefix,"Nsig","svg",sep="."),width=5,height=4)
	hist(nsig,breaks=40,xlab="Differentially expressed genes (FDR<0.05)",ylab="Number of comparisons",col="grey",main=titlePrefix)
	abline(v=mean(nsig),lty="dashed",lwd=2,col="yellow")                                         
	legend("topright",legend=c(paste("Mean:",round(mean(nsig)),"DEGs"),paste("Median:",round(median(nsig)),"DEGs")),bty="n")
        dev.off()
        svg(paste(saveToFilePrefix,"logFC","svg",sep="."),width=5,height=4)
	hist(res$mean.logFC.below.p05,breaks=40,xlab="average log2(foldchange) at FDR<0.05",ylab="Number of comparisons",col="grey",main=titlePrefix)
        dev.off()
	if(addDE){
		svg(paste(saveToFilePrefix,"sensitivity","svg",sep="."),width=5,height=5)
		getSensitivityMatrix(readPermResults(permres))
		dev.off()
	}
    }
    if(returnResults) return(permres)
}

.ttp <- function(x,groups,paired=F){
	x <- as.numeric(x)
	if(paired){
		d <- data.frame(indiv=rep(1:(length(x)/2),2),group=groups,x=x)
		a <- try(aov(x~indiv+group,data=d),silent=T)
		if("aov" %in% class(a)){
			return(summary(a)[[1]][["Pr(>F)"]][[2]])
		}
	}else{
		t <- try(t.test(x[which(groups==unique(groups)[[1]])],x[which(groups==unique(groups)[[2]])]),silent=T)
		if(class(t)!="try-error")	return(t$p.value)
	}
	return(NA)
}

.doPhenotest <- function(X, pheno, testfunc=.ttp, paired=F, progBar=NULL, progFile=NULL){
    if(is.null(X)) return(NULL)
    if(class(X)=="character" & length(X)==1){
        e <- pheno[strsplit(as.character(X),",",fixed=T)[[1]],]
    }else{
        e <- X
    }
    res <- as.numeric(apply(e,2,groups=rep(c("A","B"),each=nrow(e)/2),paired=paired,FUN=testfunc))
    if(!is.null(progFile) & !is.null(progBar)){
        prog <- as.numeric(readBin(progFile,"double"))+1
        writeBin(prog,progFile)
        setTxtProgressBar(progBar, prog)
        flush.console()
    }
    return(res)
}

#' cellpheno.permutateIndividuals
#'
#' Tests for differences in iPSC cellular morphology between groups of random individuals.
#'
#' @param nbIndividuals A positive integer indicating the number of individuals in each group (default 2).
#' @param nbClone The number of iPSC clones to consider from each individual (default 2)
#' @param maxTests The maximum number of tests to run (default 300).
#' @param logT Whether to log-transform the data before testing (default FALSE)
#' @param testfunc Function used for testing. Defaults to `.ttp` (t-test)
#' @param agTechRep Function used to aggregate technical replicates (default: mean)
#' @param doSave Whether to save the results (default TRUE).
#' @param returnResults Whether to return the results (default FALSE).
#' @param ncores Integer; number of cores to use. Defaults to detected number of cores minus one.
#'
#' @return Nothing (results save to files), or a data.frame if returnResults=TRUE.
#'
#' @export
cellpheno.permutateIndividuals <- function(nbIndividuals=2, nbClone=2, maxTests=300, logT=FALSE, testfunc=.ttp, agTechRep=mean, doSave=TRUE, returnResults=FALSE, ncores=NULL){
	if(length(nbClone) != 1 | !(nbClone %in% c(1,2)))  stop("nbClone should be either 1 or 2.")
	if(length(nbIndividuals)!=1 | !(nbIndividuals > 1)) stop("nbIndividuals should be an integer greater than 1")
	saveToFilePrefix <- paste("cellpheno.",nbIndividuals,"indiv.",nbClone,sep="")
	data("cellpheno")
	pheno <- aggregate(cellpheno[,-1:-2],by=list(clone=cellpheno$clone),FUN=agTechRep)
	row.names(pheno) <- pheno$clone
	pheno$clone <- NULL
	w <- !duplicated(cellpheno[,1])
	anno <- data.frame(row.names=cellpheno[w,1],individual=cellpheno[w,2])
	
	if(logT) pheno <- log(pheno)
	m <- combn(unique(as.character(anno$individual)),2*(nbIndividuals*nbClone))
	if(ncol(m)>maxTests)	m <- m[,sample(1:ncol(m),maxTests)]
	comb <- apply(m,2,anno=anno,nbClone=nbClone,FUN=function(x,nbClone,anno){ 
			paste(sapply(as.character(x),anno=anno,nbClone=nbClone,FUN=function(y,nbClone,anno){ 
				y <- row.names(anno)[which(anno$individual==y)]
				if(nbClone==1)	return(sample(y,1))
				return(paste(y,collapse=","))
			}), collapse=",")
		})

	res <- t(sapply(as.character(comb), pheno=pheno, testfunc=testfunc, FUN=.doPhenotest))
	colnames(res) <- colnames(pheno)

    if(doSave)	save(res,file=paste(saveToFilePrefix,"RData",sep="."))

    d <- as.data.frame(-log10(res))
    if(!.checkPkg("beanplot")){
            warning("Skipping plots: package 'beanplot' must be installed for plotting the results...")
    }else{
	library(beanplot)
	if(doSave) png(paste(saveToFilePrefix,"png",sep="."),width=1000,height=600)
	beanplot(d,ylab="-log10(p-value) across permutations",what=c(0,1,1,1), col=c("black","black","darkgrey","red"), cutmin=0,log="",
		main=paste(nbIndividuals," individuals (",nbClone," clone) per group",sep=""))
	abline(h=-log10(0.05),lty="dashed",col="red")
	if(doSave) dev.off()
    }
    if(returnResults) return(res)
}


.randomDEgenes <- function(ngenes, foldchanges, nbPerFC){
	# select [nbPerFC] ranges of genes expression (excluding the top 10 genes and the bottom ones), 
	# then in each of these ranges select 1 gene for each foldchange (and 1/foldchange)
	x <- sapply( floor((ngenes-60)/(nbPerFC+3))*1:nbPerFC+60, FUN=function(x){ sample( (x-50):(x+50), 2*length(foldchanges)) })
	# increment genes selected twice until there are no more replicates
	x <- as.numeric(x)
	while(any(duplicated(x)))	x[duplicated] <- x[duplicated]+1
	return(x)
}

.doDEA <- function(X, mm=NULL, DEgenes=NULL, DEfcs=NULL, DEsamples=NULL, paired=FALSE, nested=NULL, DEAfunc=edgeRwrapper, progBar=NULL, progFile=NULL){
    if(is.null(X)) return(NULL)
    if(class(X)=="character" & length(X)==1){
        e <- agc[,strsplit(as.character(X),",",fixed=T)[[1]]]
    }else{
        e <- X
    }
    if(!is.null(nested)) nested <- nested[colnames(e)]
    if(!is.null(DEgenes) & !is.null(DEfcs)){
        for(i in DEsamples) e[DEgenes,i] <- e[DEgenes,i]*DEfcs
    }
    res <- DEAfunc(e, mm, DEsamples=DEsamples, paired=paired, nested=nested)
    if(!is.null(progFile) & !is.null(progBar)){
        prog <- as.numeric(readBin(progFile,"double"))+1
        writeBin(prog,progFile)
        setTxtProgressBar(progBar, prog)
        flush.console()
    }
    return(res)
}

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




#' readPermResults
#'
#' @param resFiles character vector indicating the names of the results files (.RData files save by the DEA permutation functions)
#' @param threshold the FDR threshold below which a test is considered positive (default 0.05).
#'
#' @return A list with the following slots:
#' `nbComps: number of comparisons.
#' `FP: number of false (or spurious) positives for each comparison.
#' `TP: (only if the permutations were generated with addDE=TRUE) The number of true positives for each comparison.
#' `DEGs: (only if the permutations were generated with addDE=TRUE) a data.frame containing, for each gene, the number of times it was found with a FDR below the threshold, the absLog2FC introduced in the input, and the mean fragment count.
#'
#' @export
readPermResults <- function(resFiles, threshold=0.05){
	if(class(resFiles)=="character"){
	  ll <- lapply(resFiles, threshold=threshold, FUN=function(x,threshold){
	    o <- load(x)
	    readPermResults(get(o[[1]]),threshold)
	  })
	  names(ll) <- gsub(".RData","",basename(resFiles),fixed=T)
	  return(ll)
	}else{
	  permres <- resFiles
          la <- list(nbComps=ncol(permres$PValues))
          if("inputDE.foldchange" %in% colnames(permres$byGene)){
            g <- row.names(permres$byGene)[which(permres$byGene$inputDE.foldchange==1)]
            la[["FP"]] <- apply(permres$FDR[g,],2,threshold=threshold,FUN=function(x,threshold){ sum(x < threshold) })
            g <- row.names(permres$byGene)[which(permres$byGene$inputDE.foldchange!=1)]
            la[["TP"]] <- apply(permres$FDR[g,],2,threshold=threshold,FUN=function(x,threshold){ sum(x < threshold) })
            la[["DEGs"]] <- permres$byGene[g,c("FDR.below.05","inputDE.foldchange","meanCount")]
            if(threshold != 0.05){
		la[["DEGs"]][,1] <- apply(permres$FDR[g,],1,threshold=threshold,FUN=function(x,threshold){ sum(x < threshold) })
            }
            la[["DEGs"]][,2] <- abs(log2(la[["DEGs"]][,2]))
            la[["DEGs"]][,3] <- log(la[["DEGs"]][,3])
            names(la[["DEGs"]]) <- c("FDR.below.threshold","absLog2FC","logMeanCount")
          }else{
            la[["FP"]] <- apply(permres$FDR,2,threshold=threshold,FUN=function(x,threshold){ sum(x < threshold) })
          }
          return(la)
	}
}

.getSexRatio <- function(x,id){
    tt <- table(as.character(id[x]))
    if(length(tt)==1)   return(names(tt))
    tt[[1]]/tt[[2]]
}


#' getSensitivityMatrix
#'
#' Retrieves and eventually plot a sensitivity matrix across foldchanges and expression levels
#'
#' @param res A list such as one of the items returned by the `readPermResults()` function.
#' @param bins Either an integer indicating the number of bins in which to (try to) split the expression levels (default 5), or a list of ranges for binning.
#' @param unlogExpr Whether to "un-log" read counts for labels (default TRUE).
#' @param doPlot Whether to plot rather than return the matrix (default TRUE).
#'
#' @return either a heatmap or a matrix.
#'
#' @export
getSensitivityMatrix <- function(res, bins=5, unlogExpr=T, doPlot=T){
    if(!.checkPkg("pheatmap") | !.checkPkg("grid"))	stop("The 'grid' and 'pheatmap' packages must be installed in order to plot the sensitivity matrix.")
    d <- res$DEGs[order(res$DEGs$logMeanCount),]
    d$FC <- as.character(2^d$absLog2FC)
    if(length(bins)==1){
	cc <- as.character(cut(d$logMeanCount, bins))
    }else{
	cc <- sapply(d$logMeanCount,bins=bins,FUN=function(g,bins){
		which(sapply(bins,g=g,FUN=function(x,g){
			g >= x[[1]] & g <= x[[2]]
		}))[[1]]
	})
	cc <- sapply(bins,collapse=",",FUN=paste)[cc]
    }
    fcs <- sort(unique(d$FC),decreasing=T)
    m <- matrix(0, nrow=length(fcs),ncol=length(unique(cc)))
    colnames(m) <- unique(cc)
    rownames(m) <- fcs
    for(i in 1:nrow(m)){
        for(j in 1:ncol(m)){
            w <- which(cc==colnames(m)[j] & d$FC==rownames(m)[i])
            m[i,j] <- sum(d$FDR.below.threshold[w])/(length(w)*res$nbComps)
        }
    }
    if(unlogExpr)	colnames(m) <- sapply(colnames(m),FUN=function(x){ paste(round(exp(as.numeric(strsplit(gsub("]","",gsub("(","",x,fixed=T),fixed=T),",")[[1]]))-1),collapse="-") })
    if(doPlot){
        require(pheatmap)
        require(grid)
        setHook("grid.newpage", function() pushViewport(viewport(x=0.95,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
        pheatmap(100*m,col=colorRampPalette(c("white","red","darkred"))(100),cluster_rows=F,cluster_cols=F,display_numbers=T,number_format="%.0f %%",number_color="black",border_color=NA,main="Sensitivity by DEGs' foldchange and expression level")
	setHook("grid.newpage", NULL, "replace")
	grid.text(ifelse(unlogExpr,"Expression range (median read count)","Expression range, in log(read count)"), y=-0.07, gp=gpar(fontsize=16))
	grid.text("Foldchange", x=1, rot=90, gp=gpar(fontsize=16))
    }else{
        return(m)
    }
}

#' getSensitivityMatrices
#'
#' Plot a joined sensitivity matrix across foldchanges and expression levels
#'
#' @param reslist A list of analysis results, such as produced by the `readPermResults()` function.
#' @param bins Either an integer indicating the number of bins in which to (try to) split the expression levels (default 5), or a list of ranges for binning.
#' @param unlogExpr Whether to "un-log" read counts for labels (default TRUE).
#'
#' @return a heatmap.
#'
#' @export
getSensitivityMatrices <- function(reslist, bins=5, unlogExpr=T, display_numbers=T){
    if(!.checkPkg("pheatmap"))	stop("The 'pheatmap' package must be installed in order to plot the sensitivity matrix.")
    ml <- lapply(reslist,doPlot=F, bins=bins, unlogExpr=unlogExpr,FUN=getSensitivityMatrix)
    m <- ml[[1]]
    if(length(ml)>1){
        for(i in 2:length(ml)){
            m <- cbind(m,ml[[i]])
        }
    }
    library(pheatmap)
    pheatmap(100*m,col=colorRampPalette(c("white","red","darkred"))(100),cluster_rows=F,cluster_cols=F,display_numbers=display_numbers,number_format="%.0f%%",number_color="black",border_color=NA,gaps_col=cumsum(lapply(ml,FUN=ncol)))
}

#' TMM
#'
#' A wrapper for edgeR's TMM normalization
#'
#' @param dat The data to be normalized.
#'
#' @return A dataframe of normalized read counts.
#'
#' @export
TMM <- function(dat){
    library(edgeR)
    d <- DGEList(counts = dat, group = colnames(dat))
    d <- calcNormFactors(d)
    return(as.data.frame(cpm(d, normalized.lib.sizes = T) * median(d$samples$lib.size)/1e+06))
}



#' getRecurrentGenes
#'
#' Returns the frequencing of each gene being found differentially expressed in the permutation analyses.
#'
#' @param permResults Either the results of the `DEA.permutateIndividuals` or `DEA.permutateClones` function, or a character vector of the files in which the outpout of these functions was save.
#' @param removeInputDE Whether to remove the input DEGs (default TRUE)
#' @param varToFetch Variable to fetch, any of: `pval.below.05`, `pval.below.01`, `FDR.below.05`, `FDR.below.01` (default `pval.below.05`).
#' @param d For internal use.
#'
#' @return A data.frame.
#'
#' @export
getRecurrentGenes <- function(permResults, removeInputDE=TRUE, varToFetch="pval.below.05", d=NULL){
	if(class(permResults)=="character"){
		i <- 0
		if(!is.null(d)) i <- ncol(d)
		for(f in permResults){
			load(f)
			d <- getRecurrentGenes(permres, removeInputDE=removeInputDE, varToFetch=varToFetch, d=d)
		}
		colnames(d)[(i+1):ncol(d)] <- gsub(".RData","",permResults,fixed=T)
		return(d)
	}
	bg <- permResults$byGene
	if(removeInputDE & "inputDE.foldchange" %in% names(bg))	bg[which(bg$inputDE.foldchange!=1),varToFetch] <- NA
	if(is.null(d))	d <- data.frame(row.names=row.names(bg))
	cname <- paste("V",ncol(d)+1,sep="")
	d[[cname]] <- NA
	d[row.names(bg),cname] <- bg[[varToFetch]]
	return(d)
}


#' getSpuriousFCs
#'
#' Returns the log2(foldchange) of all spurious DEGs in results of permutations analyses.
#'
#' @param permResults Either the results of the `DEA.permutateIndividuals` or `DEA.permutateClones` function, or a character vector of the files in which the output of these functions was saved.
#' @param pthreshold The p-value threshold above which to consider a gene differentially-expressed.
#'
#' @return A vector of log2 foldchanges.
#'
#' @export
getSpuriousFCs <- function(permResults, pthreshold=0.05){
	if(class(permResults)=="character"){
		res <- c()
		for(f in permResults){
			load(f)
			res <- c(res, getSpuriousFCs(permres,pthreshold=pthreshold))
		}
		return(res)
	}
	fcs <- permResults$logFC
	pv <- permResults$PValues
	if("inputDE.foldchange" %in% names(permResults$byGene)){
		w <- row.names(permResults$byGene)[permResults$byGene$inputDE.foldchange==1]
		pv <- pv[w,]
		fcs <- fcs[w,]
	}
	return(round(as.numeric(as.matrix(fcs)[which(as.matrix(pv)<pthreshold)]),2))
}


#' getPermROC
#'
#' Plots the ROC curve of permutation DEAs results
#'
#' @param permResults Either the results of the `DEA.permutateIndividuals` or `DEA.permutateClones` function, or a character vector of length 1 indicating the files in which such an object was saved.
#' @param minPts Integer>0 indicating the minimum approximative number of points over which to plot the curve (default 100). Higher values produce a higher resolution plot, but will take longer to run.
#' @param qprobs Numeric vector of length 2 indicating the probabilities for the confidence interval (default `c(0.05,0.95)`). Set to NA to disable.
#' @param main Plot title (default "ROC Curve")
#' @param col Color (default blue)
#' @param type Plot type (default 'both')
#' @param lty Line type (default 1)
#' @param add Logical; whether to add the curve on an existing plot (default FALSE)
#' @param showAUC Logical; whether to display the area under the curve (default TRUE, unless add=TRUE)
#' @param returnCoords Logical; whether to return the point coordinates instead of plotting them (default FALSE)
#' @param ncores Integer; number of cores to use (default 1)
#'
#' @return Produces a plot, or returns a list with the slots `sens` and `spec` if returnCoords=TRUE.
#'
#' @export
getPermROC <- function(permResults, minPts=100, qprobs=c(0.05,0.95), main=NULL, col="blue", type="b", lty=1, add=F, lwd=2, showAUC=!add, returnCoords=FALSE, pch=1, ncores=1){
	if(class(permResults)=="character"){
		load(permResults)
		permResults <- permres
		rm(permres)
	}
	atP <- sort(.getMinPts(as.numeric(as.matrix(permResults$PValues)),minPts=minPts))
	sig <- permResults$byGene[row.names(permResults$PValues),"inputDE.foldchange"] != 1
	
	cl <- initPar(ncores)
	if(any(class(cl) != "logical")){
		clusterExport(cl, c(".sens",".spec"), envir=environment())
		sens <- t(parSapply(cl,atP,sig=sig,p=permResults$PValues,qprobs=qprobs, FUN=function(x, sig, p, qprobs){ 
			tmp <- apply(p, 2, sig=sig, thres=x, FUN=.sens)
			return(quantile(tmp,c(qprobs[1],0.5,qprobs[2])))
		}))
		spec <- t(parSapply(cl,atP,sig=sig,p=permResults$PValues,qprobs=qprobs, FUN=function(x, sig, p, qprobs){ 
			tmp <- apply(p, 2, sig=sig, thres=x, FUN=.spec)
			return(quantile(tmp,c(qprobs[1],0.5,qprobs[2])))
		}))
		stopCluster(cl)
	}else{
		sens <- t(sapply(atP,sig=sig,p=permResults$PValues,qprobs=qprobs, FUN=function(x, sig, p, qprobs){ 
			tmp <- apply(p, 2, sig=sig, thres=x, FUN=.sens)
			return(quantile(tmp,c(qprobs[1],0.5,qprobs[2])))
		}))
		spec <- t(sapply(atP,sig=sig,p=permResults$PValues,qprobs=qprobs, FUN=function(x, sig, p, qprobs){ 
			tmp <- apply(p, 2, sig=sig, thres=x, FUN=.spec)
			return(quantile(tmp,c(qprobs[1],0.5,qprobs[2])))
		}))
	}

	w <- which(!duplicated(cbind(sens[,2],spec[,2])))
	sens <- sens[w,]
	spec <- spec[w,]
	row.names(sens) <- atP
	row.names(spec) <- atP
	
	if(returnCoords) return(list(sens=sens, spec=spec))
	
	if(add){
		lines(1-spec[,2], sens[,2],type=type,lty=lty,lwd=lwd,col=col,pch=pch)
	}else{
		plot(1-spec[,2], sens[,2],xlab="1-Specificity",ylab="Sensitivity",xlim=c(0,1),ylim=c(0,1),type=type,lty=lty,lwd=lwd,main=ifelse(is.null(main),"ROC curve",main),col=col,pch=pch)
	}
	if(class(qprobs)=="numeric" & length(qprobs)==2){
		tcol <- col2rgb(col)
		tcol <- rgb(tcol["red", 1][[1]], tcol["green", 1][[1]], tcol["blue", 1][[1]], 50, maxColorValue = 255)
		polygon(c(1-spec[,1],rev(1-spec[,3])),c(sens[,1],rev(sens[,3])),border=NA,col=tcol)
	}
	if(showAUC) legend("bottomright",bty="n",legend=c(paste("AUC:",round(auc(1-spec[,2],sens[,2]),3))))
}

.sens <- function(x, sig, thres){
	sum(sig & x <= thres)/sum(sig)
}
.spec <-  function(x, sig, thres){
	sum(!(x <= thres | sig))/sum(!sig)
}

.getMinPts <- function(pv, minPts=1000){
	if(minPts < 3 | round(minPts) != minPts)	stop("Invalid value for minPts; should be an integer > 0")
	x <- c()
	i <- 0
	while(length(x) < minPts){
		i <- i+1
		x <- unique(2^round(log2(pv+1),i)-1)
	}
	return(x)
}

#' auc
#'
#' Computes the area under a curve.
#'
#' @param x x values.
#' @param y corresponding y values.
#' @param dens Number of points.
#'
#' @return The area under the curve.
#'
#' @export
auc <- function (x, y, dens = 100){
    if(dens > length(x)) dens <- length(x)
    w <- which(!is.na(x) & !is.na(y))
    y <- y[w]
    x <- x[w]
    o <- order(x)
    x <- x[o]
    y <- y[o]
    idx = 2:length(x)
    x <- as.vector(apply(cbind(x[idx - 1], x[idx]), 1, function(x) seq(x[1], x[2], length.out = dens)))
    y <- as.vector(apply(cbind(y[idx - 1], y[idx]), 1, function(x) seq(x[1], x[2], length.out = dens)))
    idx = 2:length(x)
    as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
}


#' combinePermResults
#'
#' Merges the respective elements of two lists of permutation analyses (i.e. output of readPermResults).
#'
#' @param list1 The first list of permutation results.
#' @param list2 The second list of permutation results.
#' @param balance Logical; whether to subsample the largest of the two elements being merged so that the merge is balanced (default TRUE).
#'
#' @return A list of merged permutation results.
#'
#' @export
combinePermResults <- function(list1, list2, balance=TRUE){
	ll <- list()
	for(i in 1:min(length(list1),length(list2))){
		ll[[i]] <- .combinePermRes(list1[[i]],list2[[i]])
	}
	return(ll)
}

.combinePermRes <- function(p1,p2,balance=TRUE){
	if(balance & p1$nbComps!=p2$nbComps){
		n <- min(p1$nbComps, p2$nbComps)
		p1$FP <- sample(p1$FP, n, replace=F)
		p1$TP <- sample(p1$TP, n, replace=F)
		p2$FP <- sample(p2$FP, n, replace=F)
		p2$TP <- sample(p2$TP, n, replace=F)
		p1$DEGs[,1] <- round(n/p1$nbComps * p1$DEGs[,1])
		p2$DEGs[,1] <- round(n/p2$nbComps * p2$DEGs[,1])
	}else{
		n <- (p1$nbComps+p2$nbComps)/2
	}
	ll <- list(nbComps=2*n, FP=c(as.numeric(p1$FP),as.numeric(p2$FP)), TP=c(as.numeric(p1$TP),as.numeric(p2$TP)))
	p1$DEGs[,1] <- round((p1$DEGs[,1]/p1$nbComps)*ll$nbComps)
	p2$DEGs[,1] <- round((p2$DEGs[,1]/p2$nbComps)*ll$nbComps)
	ll$DEGs <- rbind(p1$DEGs,p2$DEGs)
	return(ll)
}


#' plotPrecisionScatter
#'
#' Plots the sensitivity and false discovery rate of a list of results of permutation analyses (i.e. permres files or output of readPermResults).
#'
#' @param results A vector of permutation results files or the output of `readPermResults`.
#' @param plotSD Logical; whether to plot the standard deviation for each point (default TRUE).
#' @param x X axis value, either "fp" (false positives) or "fdr" (false discovery rate, default).
#' @param labels An optional character vector of labels to print.
#' @param acol color of the SD bars (default gray)
#' @param alty line type of the SD bars (default solid)
#' @param alwd line width of the SD bars (default 1)
#' @param ... arguments passed to the plot function.
#'
#' @return Nothing, but generates a plot.
#'
#' @export
plotPrecisionScatter <- function(results, plotSD=T, x="fdr", labels=NULL, acol="grey", alty=1, alwd=1, threshold=0.05, ...){
	x <- match.arg(tolower(x),c("fp","fdr"))
	if(class(results)=="character")	results <- readPermResults(results, threshold=threshold)
	msens <- sapply(results, FUN=function(x){ mean(x$TP)/nrow(x$DEGs) })
	if(x=="fdr"){
		mfpr <- sapply(results, FUN=function(x){ mean(x$FP/(x$FP+x$TP)) })
	}else{
		mfpr <- sapply(results, FUN=function(x){ mean(x$FP) })
	}	
	
	xlab <- ifelse(x=="fdr","False/spurious discovery rate","False/spurious DEGs")
	if(plotSD){
		sensd <- sapply(results, FUN=function(x){ sd(x$TP/nrow(x$DEGs)) })
		if(x=="fdr"){
			fprd <- sapply(results, FUN=function(x){ sd(x$FP/(x$FP+x$TP)) })
		}else{
			fprd <- sapply(results, FUN=function(x){ sd(x$FP) })
		}			
		plot(mfpr, msens, xlim=range(c(mfpr-fprd, mfpr+fprd)), ylim=range(c(msens-sensd,msens+sensd)), ylab="Sensitivity", xlab=xlab, ...)
		arrows(mfpr, msens-sensd, mfpr, msens+sensd, length=0.05, angle=90, code=3, col=acol, lty=alty, lwd=alwd)
		arrows(mfpr-fprd, msens, mfpr+fprd, msens, length=0.05, angle=90, code=3, col=acol, lty=alty, lwd=alwd)
	}else{
		plot(mfpr, msens, ylab="Sensitivity", xlab=xlab, ...)
	}
	if(!is.null(labels)){
		text(mfpr, msens, labels=labels)
	}
}