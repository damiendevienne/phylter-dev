# Independent verification using ape's Newick parser and tree distances.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1,2))
root <- normalizePath(args[1])
expanded.root <- if(length(args)==2) normalizePath(args[2]) else file.path(root,"carnivora-tip-scaling")
datasets <- read.delim(file.path(expanded.root,"datasets.tsv"))
library(ape)
max.error <- 0
for (genes in unique(datasets$genes)) {
  original <- read.tree(file.path(root,"carnivora-full",paste0("sample-",genes,".nwk")))
  for (j in which(datasets$genes==genes)) {
    factor <- datasets$factor[j]
    expanded <- read.tree(file.path(expanded.root,datasets$file[j]))
    stopifnot(length(expanded)==length(original))
    for (i in seq_along(original)) {
      labels <- original[[i]]$tip.label
      stopifnot(length(expanded[[i]]$tip.label)==factor*length(labels),
                all(labels %in% expanded[[i]]$tip.label),
                all(is.finite(expanded[[i]]$edge.length)),all(expanded[[i]]$edge.length>=0))
      expected <- cophenetic.phylo(original[[i]])
      actual <- cophenetic.phylo(expanded[[i]])[labels,labels,drop=FALSE]
      error <- max(abs(actual-expected))
      stopifnot(error < 1e-12 * max(1,max(abs(expected))))
      max.error <- max(max.error,error)
    }
    cat("PASS",genes,"genes,",53*factor,"taxa\n")
  }
}
cat("Maximum original-distance change:",format(max.error,digits=17),"\n")
