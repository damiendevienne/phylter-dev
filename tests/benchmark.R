# Developer benchmark against whichever R library R_LIBS selects.
# Rscript tests/benchmark.R build/phylter
args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
binary <- normalizePath(args[1])
library(phylter)
data(carnivora)
folder <- file.path(dirname(binary),"benchmark")
dir.create(folder,showWarnings=FALSE)
input <- file.path(folder,"carnivora.nwk")
ape::write.tree(carnivora,input,digits=17)
reference <- numeric(5)
native <- numeric(5)
for (i in seq_len(5)) {
  gc(); set.seed(42)
  reference[i] <- system.time(invisible(capture.output(
    result <- phylter(carnivora,parallel=FALSE,verbose=FALSE))))[["elapsed"]]
  native[i] <- system.time(status <- system2(binary,
    c("run","--trees",shQuote(input),"--out",shQuote(file.path(folder,"native")),"--force"),
    stdout=file.path(folder,"native.stdout"),stderr=file.path(folder,"native.stderr")))[["elapsed"]]
  stopifnot(status==0)
  scores <- read.delim(file.path(folder,"native.scores.tsv"))$quality
  stopifnot(max(abs(scores-result$Final$AllOptiScores))<1e-10)
}
cat("R library:",find.package("phylter"),"\n")
cat("R seconds:",reference,"; median:",median(reference),"\n")
cat("C++ seconds (including input/output):",native,"; median:",median(native),"\n")
cat("Median ratio R/C++:",median(reference)/median(native),"\n")
