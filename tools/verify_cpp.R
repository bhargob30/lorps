#!/usr/bin/env Rscript
# Run from the repository root before a release.
expected <- strsplit(readLines("inst/CORE_SOURCE_MD5", warn = FALSE)[1L], " +")[[1L]][1L]
actual <- unname(tools::md5sum("src/finite_stickbreak_eff.cpp"))
stopifnot(identical(actual, expected))
cat("Core C++ checksum verified.\n")
