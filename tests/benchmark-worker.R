# Isolated original-R process used by benchmark-resources.py.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
library(phylter, lib.loc = args[1])
options(digits = 17)
input <- args[2]
prefix <- args[3]
trees <- ape::read.tree(input)
if (inherits(trees, "phylo")) trees <- list(trees)
names(trees) <- paste0(tools::file_path_sans_ext(basename(input)), ":", seq_along(trees))
set.seed(42)
elapsed <- system.time(invisible(capture.output(
  result <- phylter(trees, parallel = FALSE, verbose = FALSE)
)))[["elapsed"]]
outliers <- result$Final$Outliers
if (is.null(outliers)) outliers <- matrix(character(), 0, 2)
colnames(outliers) <- c("gene", "species")
write.table(outliers, paste0(prefix, ".outliers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
discarded <- result$Final$Discarded
if (is.null(discarded)) discarded <- matrix(character(), 0, 2)
colnames(discarded) <- c("gene", "species")
write.table(discarded, paste0(prefix, ".discarded.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(state = seq_along(result$Final$AllOptiScores) - 1,
                       quality = result$Final$AllOptiScores),
            paste0(prefix, ".scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(format(elapsed, digits = 17), paste0(prefix, ".analysis-seconds.txt"))
writeLines(c(paste("reference_library", find.package("phylter")),
             paste("version", packageVersion("phylter"))),
           paste0(prefix, ".reference.txt"))
