# Developer-only regression test. R is not needed by the executable or CTest.
# From the parent R repository:
# TMPDIR="$PWD/.audit/tmp" R_LIBS="$PWD/.audit/reference-lib:$PWD/.audit/library" \
#   Rscript tests/compare-reference.R build/phylter build/reference-validation
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
binary <- normalizePath(args[1])
root <- args[2]
dir.create(root, recursive=TRUE, showWarnings=FALSE)
root <- normalizePath(root)
library(phylter)
data(carnivora)
options(digits=17)
stopifnot(packageVersion("phylter") == "0.9.12")
cat("Reference library:", find.package("phylter"), "\n")
trees <- carnivora
names(trees) <- sprintf("g%03d", seq_along(trees))
matrices <- trees2matrices(trees)
set.seed(20260920)
missing <- lapply(trees, function(tr) ape::drop.tip(tr, sample(tr$tip.label,5)))
zero <- matrices[[1]] * 0
set.seed(71)
indefinite <- lapply(1:7,function(i) {
  x <- matrix(runif(64),8,8); x <- x+t(x); diag(x) <- 0
  dimnames(x) <- list(letters[1:8],letters[1:8]); x
})
names(indefinite) <- sprintf("g%03d",seq_along(indefinite))
never_together <- lapply(seq_along(matrices),function(i) {
  keep <- setdiff(seq_len(53),if(i <= 60) 1 else 2)
  matrices[[i]][keep,keep]
})
names(never_together) <- names(matrices)
supported <- trees
for (i in seq_along(supported)) supported[[i]]$node.label <- as.character(
  rep(c(95,60,100,80),length.out=supported[[i]]$Nnode))
cases <- list(
  default=list(X=matrices),
  trees=list(X=trees),
  nodal=list(X=trees, distance="nodal"),
  mean=list(X=matrices, Norm="mean"),
  unscaled=list(X=matrices, Norm="none"),
  no_islands=list(X=matrices, test.island=FALSE),
  column_normalization=list(X=matrices, normalizeby="col"),
  no_wr_normalization=list(X=matrices, normalizeby="none"),
  gene_threshold=list(X=matrices, k2=1.5),
  missing_taxa=list(X=trees2matrices(missing)),
  discarded=list(X=c(list(g000=zero),matrices)),
  initial=list(X=matrices, InitialOnly=TRUE),
  no_outliers=list(X=matrices, k=1e10),
  rejected=list(X=matrices, stop.criteria=1),
  indefinite=list(X=indefinite),
  never_together=list(X=never_together),
  support=list(X=supported,bvalue=70),
  support_nodal=list(X=supported,bvalue=70,distance="nodal")
)
trace("DistatisFast", where=asNamespace("phylter"), print=FALSE,
      exit=quote(.GlobalEnv$reference_states[[length(.GlobalEnv$reference_states)+1L]] <- returnValue()))
max_error <- 0
numeric_equal <- function(actual, expected, label) {
  stopifnot(identical(dim(actual),dim(expected)), length(actual)==length(expected),
            all(is.finite(actual)), all(is.finite(expected)))
  error <- max(abs(actual-expected),0)
  max_error <<- max(max_error,error)
  if (error > 1e-8 * max(1,max(abs(expected)))) stop(label, ": error = ", error)
}
read_numeric <- function(path) as.matrix(read.table(path,header=FALSE,check.names=FALSE))
pair_equal <- function(path,expected) {
  actual <- as.matrix(read.delim(path,check.names=FALSE,colClasses="character"))
  storage.mode(actual) <- "character"
  if (is.null(expected)) expected <- matrix(character(),0,2)
  stopifnot(identical(unname(actual),unname(expected)))
}
for (name in names(cases)) {
  case <- cases[[name]]
  folder <- file.path(root,name)
  if (dir.exists(folder)) stop("Use a fresh validation directory: ", folder)
  dir.create(folder)
  input <- file.path(folder,"input")
  dir.create(input)
  tree_input <- inherits(case$X[[1]],"phylo")
  for (gene in names(case$X)) {
    if (tree_input) ape::write.tree(case$X[[gene]],file=file.path(input,paste0(gene,".nwk")),digits=17)
    else {
      out <- data.frame(taxon=rownames(case$X[[gene]]),case$X[[gene]],check.names=FALSE)
      write.table(out,file.path(input,paste0(gene,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
    }
  }
  # Match the CLI's sorted file order, also for discarded genes.
  case$X <- case$X[order(names(case$X))]
  reference_states <- list()
  set.seed(42)
  timing_r <- system.time(ref <- do.call(phylter,modifyList(list(parallel=FALSE,verbose=FALSE),case)))[["elapsed"]]
  cli <- c("run",if(tree_input) "--trees" else "--matrices",shQuote(input),
           "--out",shQuote(file.path(folder,"result")),
           "--diagnostics",shQuote(file.path(folder,"states")))
  mapping <- c(distance="--distance",Norm="--norm",normalizeby="--normalize-by",k2="--k2",k="--k",stop.criteria="--stop",bvalue="--support-cutoff")
  for (option in names(mapping)) if (!is.null(case[[option]])) cli <- c(cli,mapping[[option]],as.character(case[[option]]))
  if (identical(case$test.island,FALSE)) cli <- c(cli,"--no-islands")
  if (isTRUE(case$InitialOnly)) cli <- c(cli,"--initial-only")
  timing_cpp <- system.time(status <- system2(binary,cli,
       stdout=file.path(folder,"stdout.log"),stderr=file.path(folder,"stderr.log")))[["elapsed"]]
  stopifnot(status==0)
  if (isTRUE(case$InitialOnly)) {
    scores <- reference_states[[1]]$quality
    expected_outliers <- expected_discarded <- expected_cells <- NULL
    species_order <- rownames(ref$WR)
    gene_order <- colnames(ref$WR)
  } else {
    scores <- ref$Final$AllOptiScores
    expected_outliers <- ref$Final$Outliers
    expected_discarded <- ref$Final$Discarded
    expected_cells <- ref$Final$CELLSREMOVED
    species_order <- ref$Final$species.order
    gene_order <- colnames(ref$Final$WR)
  }
  numeric_equal(read.delim(file.path(folder,"result.scores.tsv"))$quality,scores,paste(name,"scores"))
  pair_equal(file.path(folder,"result.outliers.tsv"),expected_outliers)
  pair_equal(file.path(folder,"result.discarded.tsv"),expected_discarded)
  cell_path <- file.path(folder,"states","cells.tsv")
  actual_cells <- if(file.info(cell_path)$size) read_numeric(cell_path) else matrix(integer(),0,2)
  if (is.null(expected_cells)) expected_cells <- matrix(integer(),0,2)
  stopifnot(identical(unname(actual_cells),unname(expected_cells)),
            identical(readLines(file.path(folder,"states","taxa.txt")),species_order),
            identical(readLines(file.path(folder,"states","genes.txt")),gene_order))
  # Compare every accepted state, not just final outlier counts.
  qualities <- vapply(reference_states,function(s) s$quality,numeric(1))
  for (i in seq_along(scores)) {
    state <- reference_states[[which.min(abs(qualities-scores[i]))]]
    base <- file.path(folder,"states",paste0("state-",i-1))
    for (field in c("rv","compromise","fgram","wr","weights")) {
      expected <- switch(field,rv=state$RVmat,compromise=state$compromise,
                         fgram=tcrossprod(state$F),wr=Dist2WR(state),weights=matrix(state$alpha))
      numeric_equal(read_numeric(paste0(base,".",field,".tsv")),expected,paste(name,i,field))
    }
  }
  cat(sprintf("PASS %-22s states=%2d outliers=%3d R=%.3fs C++(+diagnostics)=%.3fs\n",
              name,length(scores),length(expected_outliers)/2,timing_r,timing_cpp))
}
untrace("DistatisFast",where=asNamespace("phylter"))
cat("Maximum absolute numeric difference:",format(max_error,digits=17),"\n")
writeLines(c(paste("Reference:",find.package("phylter")),
             paste("Maximum absolute error:",format(max_error,digits=17)),
             capture.output(sessionInfo())),file.path(root,"provenance.txt"))
