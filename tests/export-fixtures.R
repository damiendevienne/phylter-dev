# Regenerate frozen CTest fixtures from the UNMODIFIED reference installation.
# Usage: Rscript tests/export-fixtures.R tests/fixtures
args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
folder <- args[1]
dir.create(folder,recursive=TRUE,showWarnings=FALSE)
library(phylter)
data(carnivora)
options(digits=17)
names(carnivora) <- paste0("carnivora:",seq_along(carnivora))
ape::write.tree(carnivora,file=file.path(folder,"carnivora.nwk"),digits=17)
set.seed(42)
ref <- phylter(carnivora,parallel=FALSE,verbose=FALSE)
write_matrix <- function(x,name) {
  if(is.null(x)) x <- matrix(numeric(),0,2)
  write.table(x,file.path(folder,name),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}
write_matrix(ref$Final$Outliers,"carnivora.outliers.tsv")
write_matrix(ref$Final$CELLSREMOVED,"carnivora.cells.tsv")
write_matrix(ref$Final$AllOptiScores,"carnivora.scores.tsv")
for (name in c("Initial","Final")) {
  state <- ref[[name]]
  for (field in c("WR","RV","weights","compromise","FGram")) {
    x <- if(field=="FGram") tcrossprod(state$F) else state[[field]]
    write_matrix(x,paste0("carnivora.",tolower(name),".",tolower(field),".tsv"))
  }
}
set.seed(972)
samples <- c(list(1,c(1,2),rep(1,101),rep(1,100),-50:50,c(rep(0,60),1:40)),
             lapply(c(3,4,20,99,100,101,2000),function(n) exp(rnorm(n))),
             lapply(1:20,function(i) sample(0:10,200,replace=TRUE)))
con <- file(file.path(folder,"medcouple.tsv"),"w")
for (x in samples) writeLines(paste(c(format(as.numeric(medcouple(x)),digits=17),
                                    format(x,digits=17)),collapse="\t"),con)
close(con)
set.seed(822)
for (i in 1:10) {
  x <- matrix(sample(0:9,144,replace=TRUE),12,12)
  x <- x+t(x); diag(x) <- 0
  write_matrix(x,paste0("cluster-",i,".tsv"))
  write_matrix(hclust(as.dist(x))$order,paste0("cluster-",i,".order.tsv"))
}
writeLines(c("Generated from phylter 0.9.12, original revision:",
             "4d74241169be3882da9d5f08da47bd40604764eb",
             paste("Reference library:",find.package("phylter")),
             "Source data: ../data/carnivora.rda, MD5 3c8348ab430dcdacef2c2c4dcc95580d",
             "Eigenvectors compared using F %*% t(F); signs are arbitrary.",
             capture.output(sessionInfo())),file.path(folder,"PROVENANCE.txt"))
