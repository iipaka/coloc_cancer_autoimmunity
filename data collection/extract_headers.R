# ==============================================================================
# Purpose: Examine column headers in GWAS summary statistics files
# ==============================================================================

library(data.table)

gwas_summarystats <- fread("Data_raw/gwas_summarystats.csv")
header_list <- list()

#For every row of GWAS study in the csv:
for (i in 1:nrow(gwas_summarystats)) {
  s_name <- gwas_summarystats$name[i] #Extract study name
  f_path <- gwas_summarystats$file_path[i] #Extract file path
  
  all_cols <- names(fread(f_path, nrows = 0)) #Read column names from that file
  dt_row <- as.data.table(t(all_cols)) #Transpose column names (from vertical list) and convert to data table
  dt_row <- cbind(Study = s_name, dt_row) #Add new column at the start of the row displaying study name
  header_list[[i]] <- dt_row #Save row into list at position i
}

full_report <- rbindlist(header_list, fill = TRUE) #Combine all rows into master table
col_count <- ncol(full_report) - 1 #Calculate max number of columns seen across all files
new_names <- c("Study", paste0("Col_", 1:col_count)) #Create naming sequence for columns

setnames(full_report, names(full_report), new_names) #Apply new column names to master table

fwrite(full_report, "Output/full_header_report.csv") #Save as new CSV file

message("Saved header summary to: Output/full_header_report.csv")

