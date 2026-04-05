outfile <- commandArgs(trailingOnly = TRUE)[1]
writeLines(capture.output(sessionInfo()), con = outfile)
