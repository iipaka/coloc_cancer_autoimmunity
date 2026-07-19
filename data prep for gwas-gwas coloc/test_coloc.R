# ==============================================================================
# Purpose: Test ALL processed .rds files for coloc compatibility
# ==============================================================================
library(coloc)

all_files <- list.files("Data_processed/", pattern = "\\.rds$", full.names = TRUE)
message(paste("Found", length(all_files), "files to test...\n"))

bad_files <- c()

for(f in all_files) {
  my_data <- readRDS(f)
  
  check_status <- tryCatch({
    check_dataset(my_data)
    TRUE
  }, error = function(e) {
    message(paste("FAILED:", basename(f)))
    message(paste("  -> Reason:", e$message))
    return(FALSE)
  })
  
  if(!check_status) {
    bad_files <- c(bad_files, f)
  }
}

# Print the final summary
message("\n==================================================")
if(length(bad_files) == 0) {
  message("SUCCESS: All ", length(all_files), " files.")
  message("Ready for coloc")
} else {
  message("WARNING: ", length(bad_files), " files not ready for coloc.")
}
message("==================================================")