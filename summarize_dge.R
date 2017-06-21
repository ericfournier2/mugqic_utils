#' Import a set of identically structured files.
#'
#' Given a set of files with identical row names and column names, this
#' function reads all files and concatenate the requested columns from each.
#'
#' @param file.names The files to be read.
#' @param header.columns Indices or names of row-identifying columns which should
#'   be repeated across all files. Those columns are added only once to the
#'   output, as the very first columns.
#' @param data.columns Indices or names fo the columns containing unique data in
#'   each file. The values from each file will be added to the output.
#' @param file.labels A vector of labels for the imported files. This must be of
#'   the of same length as \code{file.names}. The label is prefixed to column names
#'   in the resulting data frame.
#' @param id.col The index or name of the column containing row names.
#' @return A \code{data-frame} with the concatenated information from all files.
#' @export
summarize.files <- function(file.names, header.columns, data.columns, file.labels=basename(file.names), id.col=0, ...) {
    # Read all files and get all IDs
    raw.files = list()
    id.list = c()
    for(i in 1:length(file.names)) {
        this.data = read.table(file.names[i], stringsAsFactors=FALSE, ...)
        rownames(this.data) = this.data[,id.col]
        raw.files[[file.labels[i]]] = this.data
    }
    
    # Merge recursively.
    results=raw.files[[1]][, header.columns]
    for(i in 2:length(raw.files)) {
        results = merge(results, raw.files[[i]][, header.columns], all=TRUE)
        rownames(results) = results[,id.col]
    }
    
    # Add data columns
    for(i in 1:length(raw.files)) {
        old.names = colnames(results)
        id.matches = match(rownames(results), raw.files[[i]][,id.col])
        results = cbind(results, raw.files[[i]][id.matches, data.columns])
        colnames(results) = c(old.names, paste(colnames(raw.files[[i]][,data.columns]), file.labels[i], sep="."))
    }
    
    return(results)
}

summarize.dge <- function(input.dir, dge.library="both", subdir="DGE") {
    if(dge.library=="DESeq") {
        data.columns=3:6
    } else if(dge.library=="edgeR") {
        data.columns=c(3,4,7,8)
    } else {
        data.columns=3:8
    }
    
    dge.dir = file.path(input.dir, subdir)
    dge.files = file.path(list.files(dge.dir), "dge_results.csv")
    dge.files = dge.files[file.exists(file.path(dge.dir, dge.files))]
    
    dge.names = dirname(dge.files)
    results = summarize.files(file.names=file.path(dge.dir, dge.files), 
                              header.columns=1:2, 
                              data.columns=data.columns, 
                              file.labels=dge.names, 
                              id.col="id", 
                              sep="\t", header=TRUE)
    
    return(results)
}

dge.by.comparison <- function(dge.summary, fc.threshold, p.threshold) {
    which.adj = grepl(".adj.p.value", colnames(dge.summary))
    
    results = list()
    for(i in which(which.adj)) {
        comp.name = gsub(".*.adj.p.value.(.*)", "\\1", colnames(dge.summary)[i])
        p.pass = dge.summary[,i] < p.threshold
        fc.pass = abs(dge.summary[,i - 3]) > fc.threshold
        up = dge.summary[,i - 3] > 0
        down = dge.summary[,i - 3] < 0
        
        col.select = c(2, i-2, i-3, i)
        up.select = p.pass & fc.pass & up
        gene.up = dge.summary[up.select & !is.na(up.select),col.select]
        colnames(gene.up) = gsub(paste0("(.*)\\.", comp.name), "\\1", colnames(gene.up))
        #gene.up = gene.up[!is.na(gene.up)]
        
        down.select = p.pass & fc.pass & down
        gene.down = dge.summary[down.select & !is.na(down.select),col.select]
        colnames(gene.down) = gsub(paste0("(.*)\\.", comp.name), "\\1", colnames(gene.down))
        #gene.down = gene.down[!is.na(gene.down)]
        
        results[[comp.name]] = list(Up=gene.up,
                                    Down=gene.down)
    }

    return(results)
}
