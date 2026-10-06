# Emit allowlisted UTF-8 source content, or a binary-safe path list with --paths.
args <- commandArgs(trailingOnly = TRUE)
root <- if (file.exists("DESCRIPTION")) "." else "TCGASigSurvival"
relative <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
top <- c("DESCRIPTION", "NAMESPACE", "README.md", "NEWS.md", "LICENSE", "LICENSE.md",
         "CITATION.cff", ".Rbuildignore", ".gitignore")
allowed <- relative %in% top | grepl("^(R|man|inst|tests|scripts|docs|[.]github)/", relative)
allowed <- allowed & !grepl("(^|/)([.]DS_Store|results|[.]git)(/|$)|[.](rds|Rds|log|tar[.]gz|zip)$", relative)
relative <- relative[allowed]
paths_only <- "--paths" %in% args
args <- setdiff(args, "--paths")
if (length(args)) relative <- relative[startsWith(relative, args[1])]
if (paths_only) {
  cat(jsonlite::toJSON(relative, auto_unbox = FALSE, pretty = FALSE))
  quit(status = 0)
}
binary <- grepl("[.](pdf|png|jpe?g|svgz|tiff?|gif|webp|ico|woff2?|ttf|otf)$", relative, ignore.case = TRUE)
if (any(binary)) stop("Binary documentation assets cannot be emitted as UTF-8 content. Use --paths for a complete file list and upload binaries through Git or the normal file-upload API.")
items <- lapply(relative, function(path) {
  lines <- readLines(file.path(root, path), encoding = "UTF-8", warn = FALSE)
  list(path = path, mode = "100644", type = "blob", content = paste0(paste(lines, collapse = "\n"), "\n"))
})
cat(jsonlite::toJSON(items, auto_unbox = TRUE, pretty = FALSE))
