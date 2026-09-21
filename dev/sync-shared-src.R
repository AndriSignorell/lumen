# Prüft, ob die gemeinsam genutzten C++-Header in lumen und DescToolsX
# bytegleich sind. Sie waren es nicht: lumen hatte einen älteren Stand
# von boot_framework.h ohne die Schutzabfragen.
#
# Aufruf aus einem der beiden Paketwurzelverzeichnisse, oder mit
# angepassten Pfaden.

syncCheck <- function(lumen      = "../lumen/src",
                      desctoolsx = "../DescToolsX/src",
                      files      = c("boot_framework.h", "bca_helpers.h", "extractBootArgs.R")) {

  paths <- file.path(rep(c(lumen, desctoolsx), each = length(files)), files)
  miss  <- paths[!file.exists(paths)]
  if (length(miss))
    stop("nicht gefunden: ", paste(miss, collapse = ", "))

  md5 <- unname(tools::md5sum(paths))
  dim(md5) <- c(length(files), 2L)

  out <- data.frame(file  = files,
                    lumen = md5[, 1],
                    dtx   = md5[, 2],
                    equal = md5[, 1] == md5[, 2],
                    stringsAsFactors = FALSE)

  if (!all(out$equal))
    warning("Diese Dateien laufen auseinander: ",
            paste(out$file[!out$equal], collapse = ", "))

  out
}
