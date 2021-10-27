#' getSensitivityMatrix
#'
#' Retrieves and eventually plot a sensitivity matrix across foldchanges and expression levels.
#' Use the higher-level function `getSensitivityMatrices` instead.
#'
#' @param res A list such as one of the items returned by the `readPermResults()` function.
#' @param bins Either an integer indicating the number of bins in which to (try to) split the expression levels (default 5), 
#' or a vectors of breaks defining the bins.
#' @param unlogExpr Whether to "un-log" read counts for labels (default TRUE).
#' @param doPlot Whether to plot rather than return the matrix (default TRUE).
#'
#' @return either a heatmap or a matrix.
#'
#' @export
getSensitivityMatrix <- function(res, bins=5, unlogExpr=T, doPlot=T){
  if(!.checkPkg("pheatmap") | !.checkPkg("grid"))	stop("The 'grid' and 'pheatmap' packages must be installed in order to plot the sensitivity matrix.")
  d <- res$DEGs[order(res$DEGs$logMeanCount),]
  if(unlogExpr) d$logMeanCount <- exp(d$logMeanCount)
  d$FC <- as.character(2^d$absLog2FC)
  if(length(bins)==1){
    cc <- as.character(cut(d$logMeanCount, bins))
  }else{
    cc <- as.character(cut(d$logMeanCount, unique(bins)))
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
#' Plot one or more sensitivity matrices across foldchanges and expression levels.
#' Will use `ComplexHeatmap` if available, otherwise `pheatmap`.
#'
#' @param reslist A list of analysis results, such as produced by the `readPermResults()` function.
#' The names of the element will be used as heatmap titles. Note that sensitivity can only be plotted 
#' when the analyses were run using `addDE=TRUE`.
#' @param bins Either an integer indicating the number of bins in which to (try to) 
#' split the expression levels (default 5), or a list of ranges for binning.
#' @param unlogExpr Whether to "un-log" read counts for labels (default TRUE).
#' @param do.draw Logical; whether to draw the heatmap (default TRUE). If FALSE, will return the Heatmap list object
#' (assuming ComplexHeatmap is being used - otherwise the argument is ignored)
#' @param ... Any further argument passed to the plotting function.
#'
#' @return a heatmap.
#'
#' @export
getSensitivityMatrices <- function(reslist, bins=5, unlogExpr=TRUE, display_numbers=TRUE, do.draw=TRUE, ...){
  stopifnot(is.list(reslist))
  if(all(c("nbComps", "FP") %in% names(reslist))) reslist <- list(reslist)
  if(any(sapply(reslist, FUN=function(x) is.null(x$TP))))
    stop("Some of the analyses were run without `addDE=TRUE`, and consequently sensitivity cannot be computed.")
  if(!.checkPkg("pheatmap") && !.checkPkg("ComplexHeatmap"))
    stop("One of the 'pheatmap' or 'ComplexHeatmap' package must be installed in order to plot the sensitivity matrix.")
  ml <- lapply(reslist,doPlot=F, bins=bins, unlogExpr=unlogExpr, FUN=getSensitivityMatrix)
  if(!is.null(names(ml)) && all(grepl("^clone|^[0-9]+indiv\\.vs|RData$",names(ml),ignore.case=TRUE)))
    names(ml) <- .filenameToTitle(names(ml))
  ml <- lapply(ml, FUN=function(x) round(x*100))
  if(require("ComplexHeatmap", quietly=TRUE)){
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.0f", ml[[m]][i, j]), x, y, gp = grid::gpar(fontsize = 10))
    }
    if(!display_numbers) cell_fun=NULL
    h <- NULL
    for(m in seq_along(ml)){
      h <- h + Heatmap(ml[[m]], cluster_rows=FALSE, cluster_columns=FALSE, 
                       name=ifelse(m==length(ml),"sensitivity",paste0("h",m)), cell_fun=cell_fun, 
                       show_heatmap_legend=m==length(ml), column_title=names(ml)[m], ...
              )
    }
    rt <- "Expression range (median normalized count)"
    if(length(ml)==1) rt <- gsub(" (","\n(",rt,fixed=TRUE)
    if(!do.draw) return(h)
    return(draw(h, column_title_side="bottom", column_title=rt,
         row_title_side="right", row_title="Foldchange", heatmap_legend_side="left"))
  }else{
    m <- ml[[1]]
    if(length(ml)>1){
      for(i in 2:length(ml)){
        m <- cbind(m,ml[[i]])
      }
    }
    pheatmap::pheatmap(100*m,col=colorRampPalette(c("white","red","darkred"))(100),
                       cluster_rows=F,cluster_cols=F,display_numbers=display_numbers,
                       number_format="%.0f%%",number_color="black",border_color=NA,
                       gaps_col=cumsum(lapply(ml,FUN=ncol)), ...)
  }
}


.filenameToTitle <- function(x){
  x <- gsub("\\.nested","",x)
  x <- gsub("clones\\.","",x)
  x <- gsub("Dream|Voom","",x)
  x <- gsub("\\.RData","",x)
  x <- gsub("indiv.vs.","vs",x,fixed=TRUE)
  x <- gsub("\\.1$","\n1 clone/indiv",x)
  x <- gsub("\\.2$","\n2 clones/indiv",x)
  x
}
