#!/usr/bin/env Rscript

# Run the analysis steps in order. Plotting is a separate command.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this entry point with Rscript.")
code_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
steps <- c(
  "prepare_background.R",
  "run_NC_XMR_EUR_to_AFR_CSA.R",
  "run_NC_CAMeRa_EUR_to_AFR_CSA.R",
  "run_NC_CAMeRa_EUR_AFR_FEMA_EUR.R",
  "run_NC_CAMeRa_EUR_CSA_FEMA_EUR.R",
  "run_NC_MRBEE_TL_EUR_to_AFR_CSA.R",
  "run_NC_TEMR_EUR_to_AFR_CSA.R",
  "run_NC_MRAPSS_EUR.R"
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
