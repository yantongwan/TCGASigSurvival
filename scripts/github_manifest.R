# Emit only explicitly allowlisted distribution files as UTF-8 JSON.
args <- commandArgs(trailingOnly = TRUE)
root <- if (file.exists("DESCRIPTION")) "." else "TCGASigSurvival"
relative <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
top <- c("DESCRIPTION", "NAMESPACE", "README.md", "NEWS.md", "LICENSE", "LICENSE.md",
         "CITATION.cff", ".Rbuildignore", ".gitignore")
allowed <- relative %in% top | grepl("^(R|man|inst|tests|scripts|docs|[.]github)/", relative)
allowed <- allowed & !grepl("(^|/)([.]DS_Store|results|[.]git)(/|$)|[.](rds|Rds|log|tar[.]gz|zip)$", relative)
relative <- relative[allowed]
if (length(args)) relative <- relative[startsWith(relative, args[1])]
items <- lapply(relative, function(path) {
  lines <- readLines(file.path(root, path), encoding = "UTF-8", warn = FALSE)
  list(path = path, mode = "100644", type = "blob", content = paste0(paste(lines, collapse = "\n"), "\n"))
})
cat(jsonlite::toJSON(items, auto_unbox = TRUE, pretty = FALSE))
