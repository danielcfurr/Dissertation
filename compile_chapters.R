
patterns_function <- function(file, is_appendix) {
  app <- ifelse(is_appendix, "\\\\printbibliography\\\\appendix", "")
  matrix(c(paste0("^.*\\\\input\\{", file, "-concordance}.*$"), "",
           "\\input\\{",  paste0("\\input\\{", file, "/"),
           "\\includegraphics\\{",  paste0("\\includegraphics\\{", file, "/"),
           "\\\\title\\{",          paste0(app, "\\\\newpage \\\\chapter\\{"),
           "^.*\\\\usepackage.*\\{geometry}.*$", "",
           "^.*\\\\usepackage\\{Sweave}.*$",     "",
           "^.*\\\\documentclass.*$",            "",
           "^.*\\\\bibliography\\{.*$",          "",
           "^.*\\\\author\\{.*$",                "",
           "^.*\\\\date.*$",                     "",
           "^.*\\\\begin\\{document}.*$",        "",
           "^.*\\\\end\\{document}.*$",          "",
           "^.*\\\\maketitle.*$",                "",
           "^.*\\\\printbibliography*$",         ""),
         ncol = 2, byrow = TRUE)
}

full_files  <- c("chapter_1/chapter_1.tex",
                 "chapter_2/chapter_2.tex",
                 "chapter_3/chapter_3.tex",
                 "appendix/appendix.tex")
short_files <- sub("^.+/(.+)\\.tex", "\\1", full_files)
short_files <- gsub("\\_", "\\\\_", short_files)

is_appendix <- c(F, F, F, T)

packages_lines <- c()
chapters_lines <- c()
for(j in 1:length(full_files)) {
  lines <- readLines(full_files[j])
  patterns <- patterns_function(short_files[j], is_appendix[j])
  for(i in 1:nrow(patterns))
    lines <- gsub(patterns[i,1], patterns[i,2], lines)
  chapter_line <- grep("\\chapter\\{", lines)
  packages_lines <- c(packages_lines, lines[1:(chapter_line-1)])
  chapters_lines <- c(chapters_lines, lines[chapter_line:length(lines)])
}

packages_lines <- packages_lines[packages_lines != ""]
packages_lines <- unique(packages_lines)
cat(packages_lines, sep = "\n", file = "packages.tex")

L <- length(chapters_lines)
drop <- chapters_lines[1:(L-1)] == "" & chapters_lines[2:(L)] == ""
chapters_lines <- chapters_lines[!drop]
cat(chapters_lines, sep = "\n", file = "chapters.tex")


