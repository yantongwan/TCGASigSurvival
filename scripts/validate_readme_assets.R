# Validate portable README links and byte-identical PDF assets without rendering.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else if (file.exists("DESCRIPTION")) "." else "TCGASigSurvival")
readme <- paste(readLines(file.path(root, "README.md"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
manifest <- utils::read.delim(file.path(root, "docs/figures/manifest.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 7L, !anyDuplicated(manifest$file), all(grepl("[.]pdf$", manifest$file)))
paths <- file.path(root, "docs/figures", manifest$file)
stopifnot(all(file.exists(paths)), identical(unname(tools::md5sum(paths)), manifest$md5),
  identical(as.numeric(file.info(paths)$size), as.numeric(manifest$bytes)))
for (path in paths) {
  con <- file(path, "rb"); header <- readBin(con, "raw", n = 5L); close(con)
  stopifnot(identical(rawToChar(header), "%PDF-"))
}
links <- regmatches(readme, gregexpr("\\]\\(([^)]+)\\)", readme, perl = TRUE))[[1]]
targets <- sub("^\\]\\(", "", sub("\\)$", "", links))
targets <- targets[!grepl("^(https?://|mailto:|#)", targets)]
targets <- sub("#.*$", "", targets)
stopifnot(all(file.exists(file.path(root, targets))))
for (name in manifest$file) stopifnot(grepl(paste0("docs/figures/", name), readme, fixed = TRUE))
stopifnot(!grepl("!\\[[^]]*\\]\\([^)]*[.]pdf\\)", readme, perl = TRUE))
for (path in c("README.md", "docs/figures/README.md")) {
  content <- paste(readLines(file.path(root, path), encoding = "UTF-8", warn = FALSE), collapse = "\n")
  stopifnot(!grepl("/Users/|/Volumes/|TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", content))
}
if (length(args) > 1L) {
  project <- normalizePath(args[2])
  stopifnot(identical(unname(tools::md5sum(file.path(project, manifest$source_file))), manifest$md5))
}
cat("README local links, 7 PDF headers/bytes/hashes, privacy and optional source fingerprints passed.\n")
