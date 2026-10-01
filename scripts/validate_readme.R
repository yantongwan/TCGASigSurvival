args <- commandArgs(trailingOnly = TRUE)
readme <- if (length(args)) args[1] else "README.md"
output <- if (length(args) > 1L) args[2] else file.path(tempdir(), "tcgasig_readme_validation")
lines <- readLines(readme, encoding = "UTF-8", warn = FALSE)
blocks <- list()
current <- character()
in_r <- FALSE
for (line in lines) {
  if (!in_r && identical(line, "```r")) {
    in_r <- TRUE
    current <- character()
  } else if (in_r && identical(line, "```")) {
    if ("# tutorial-demo" %in% current) blocks[[length(blocks) + 1L]] <- current
    in_r <- FALSE
  } else if (in_r) current <- c(current, line)
}
stopifnot(length(blocks) == 5L)
env <- new.env(parent = globalenv())
for (i in seq_along(blocks)) {
  cat("Running README demo block ", i, "/", length(blocks), "\n", sep = "")
  eval(parse(text = blocks[[i]]), envir = env)
  if (i == 1L) {
    env$tutorial_out <- output
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
  }
}
stopifnot(nrow(env$demo_state$results) == 2L,
          nrow(env$demo_two$continuous) == 2L,
          nrow(env$demo_two$adjusted) == 2L,
          nrow(env$demo_two$cell_high) == 2L,
          nrow(env$demo_adjusted$continuous) == 2L,
          nrow(env$demo_batch$results) == 4L,
          nrow(env$demo_protein$correlations) == 2L,
          nrow(env$demo_protein$survival) == 8L)
writeLines(c("All README tutorial checks passed.", capture.output(sessionInfo())),
           file.path(output, "README_VALIDATION.txt"))
cat("README validation completed: ", output, "\n", sep = "")
