#!/usr/bin/env Rscript

# Run the analysis steps in order. Plotting is a separate command.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this entry point with Rscript.")
code_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
steps <- c(
  "prepare_background.R",
  "run_SHBG_T2D_XMR.R",
  "run_SHBG_T2D_CAMeRa.R",
  "run_SHBG_T2D_MRAPSS.R",
  "run_SHBG_T2D_MRBEE_TL.R",
  "run_SHBG_T2D_TEMR.R"
)
args <- commandArgs(TRUE)
if (length(args)) {
  if (identical(args, "--help")) {
    cat("Analysis order:\n", paste(seq_along(steps), steps, collapse = "\n"), "\n")
    quit(status = 0)
  }
  stop("This entry point accepts only --help. Run individual scripts for method-specific options.")
}
for (step in steps) {
  message("Running ", step)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(c("--vanilla", file.path(code_dir, step))))
  if (status != 0L) stop("Analysis stopped at ", step, " (exit status ", status, ").")
}
message("Analysis complete. Run the plotting command in the analysis README.")
