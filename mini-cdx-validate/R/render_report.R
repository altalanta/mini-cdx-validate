#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript R/render_report.R <input.Rmd> [output.html]")
}

input <- args[[1]]
output <- if (length(args) >= 2) args[[2]] else {
  file.path(dirname(input), paste0(tools::file_path_sans_ext(basename(input)), ".html"))
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required")
}

render_with_rmarkdown <- function(in_file, out_file) {
  rmarkdown::render(in_file, output_file = basename(out_file), output_dir = dirname(out_file), quiet = TRUE)
}

markdown_table_to_html <- function(block) {
  trim_cells <- function(line) {
    parts <- strsplit(sub("^\\|", "", sub("\\|$", "", line)), "\\|")[[1]]
    trimws(parts)
  }
  if (length(block) < 2) {
    return(paste0("<pre>", paste(block, collapse = "\n"), "</pre>"))
  }
  header <- trim_cells(block[1])
  body <- if (length(block) > 2) block[-c(1, 2)] else character()
  rows <- lapply(body, trim_cells)
  head_html <- paste0("<thead><tr>", paste(sprintf("<th>%s</th>", header), collapse = ""), "</tr></thead>")
  body_html <- if (length(rows) > 0) {
    paste(vapply(rows, function(r) {
      paste0("<tr>", paste(sprintf("<td>%s</td>", r), collapse = ""), "</tr>")
    }, character(1)), collapse = "")
  } else {
    ""
  }
  paste0("<table>", head_html, "<tbody>", body_html, "</tbody></table>")
}

markdown_to_html <- function(lines) {
  html <- character()
  in_list <- FALSE
  i <- 1
  while (i <= length(lines)) {
    line <- lines[i]
    if (grepl("^\\s*$", line)) {
      if (in_list) {
        html <- c(html, "</ul>")
        in_list <- FALSE
      }
      i <- i + 1
      next
    }
    if (grepl("^#{1,6} ", line)) {
      level <- attr(regexpr("^#+", line), "match.length")
      text <- trimws(sub("^#+", "", line))
      html <- c(html, sprintf("<h%d>%s</h%d>", level, text, level))
      i <- i + 1
      next
    }
    if (grepl("^\\- ", line)) {
      if (!in_list) {
        html <- c(html, "<ul>")
        in_list <- TRUE
      }
      item <- trimws(sub("^\\- ", "", line))
      html <- c(html, sprintf("<li>%s</li>", item))
      i <- i + 1
      next
    }
    if (grepl("^!\\[", line)) {
      match <- regexec("^!\\[(.*)\\]\\((.*)\\)", line)
      parts <- regmatches(line, match)[[1]]
      if (length(parts) == 3) {
        alt <- parts[2]
        src <- parts[3]
      } else {
        alt <- ""
        src <- ""
      }
      fig <- sprintf("<figure><img src=\"%s\" alt=\"%s\"/><figcaption>%s</figcaption></figure>", src, alt, alt)
      html <- c(html, fig)
      i <- i + 1
      next
    }
    if (grepl("^\\|", line)) {
      block <- character()
      while (i <= length(lines) && grepl("^\\|", lines[i])) {
        block <- c(block, lines[i])
        i <- i + 1
      }
      html <- c(html, markdown_table_to_html(block))
      next
    }
    if (in_list) {
      html <- c(html, "</ul>")
      in_list <- FALSE
    }
    html <- c(html, sprintf("<p>%s</p>", line))
    i <- i + 1
  }
  if (in_list) {
    html <- c(html, "</ul>")
  }
  html
}

render_with_markdown <- function(in_file, out_file) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Fallback requires the knitr package")
  }
  tmp_md <- tempfile(fileext = ".md")
  old_root <- knitr::opts_knit$get("root.dir")
  knitr::opts_knit$set(root.dir = getwd())
  on.exit(knitr::opts_knit$set(root.dir = old_root), add = TRUE)
  knitr::knit(in_file, output = tmp_md, quiet = TRUE, envir = new.env(parent = globalenv()))
  lines <- readLines(tmp_md, warn = FALSE)
  unlink(tmp_md)
  html_body <- markdown_to_html(lines)
  cat("<html><head><meta charset=\"UTF-8\"></head><body>\n", file = out_file)
  cat(paste(html_body, collapse = "\n"), file = out_file, append = TRUE)
  cat("\n</body></html>\n", file = out_file, append = TRUE)
}

if (rmarkdown::pandoc_available("1.12.3")) {
  render_with_rmarkdown(input, output)
} else {
  message("Pandoc not found; using markdown::markdownToHTML fallback")
  render_with_markdown(input, output)
}

message("Rendered ", input, " -> ", output)
