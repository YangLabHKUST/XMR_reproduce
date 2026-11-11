#!/usr/bin/env Rscript

# Run the analysis steps in order. Plotting is a separate command.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this entry point with Rscript.")
code_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
steps <- c(
  "run_XMR_EAS_BBJ_TPMI.R",
  "run_CAMERA_EUR_to_BBJ_TPMI.R",
  "run_MRBEE_TL_EUR_to_BBJ_TPMI.R",
  "run_MRAPSS_BBJ_TPMI.R",
  "run_MRAPSS_EUR_PanUKB_5e-05.R"
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
