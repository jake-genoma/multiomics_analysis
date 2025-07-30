# helpers.R

# Helper function to read and validate data
read_and_validate <- function(file_path, header, sep, type) {
  req(file_path)
  tryCatch({
    df <- read.csv(file_path$datapath, header = header, sep = sep, 
                   stringsAsFactors = FALSE, check.names = FALSE)
    
    # Basic validation based on data type
    if (type == "transcript") {
      if (!any(grepl("Gene|ENSG", colnames(df)[1], ignore.case = TRUE))) {
        showFeedbackDanger("transcript_file", "First column should contain gene identifiers")
        return(NULL)
      }
    } else if (type == "protein") {
      if (!any(grepl("Prot|P[0-9]", colnames(df)[1], ignore.case = TRUE))) {
        showFeedbackDanger("protein_file", "First column should contain protein identifiers")
        return(NULL)
      }
    } else if (type == "metabo") {
      if (!any(grepl("Metab|HMDB|C[0-9]", colnames(df)[1], ignore.case = TRUE))) {
        showFeedbackDanger("metabo_file", "First column should contain metabolite identifiers")
        return(NULL)
      }
    } else if (type == "metadata") {
      if (colnames(df)[1] != "SampleID") {
        showFeedbackDanger("metadata_file", "First column must be 'SampleID'")
        return(NULL)
      }
    }
    
    return(df)
  }, error = function(e) {
    showNotification(paste("Error reading file:", e$message), type = "error")
    return(NULL)
  })
}

# Helper function to update metadata dropdowns
update_metadata_dropdowns <- function(session, metadata) {
  if (!is.null(metadata)) {
    meta_vars <- colnames(metadata)[-1] # exclude SampleID
    updateSelectInput(session, "corr_strata", choices = c("None", meta_vars))
    updateSelectInput(session, "diff_group", choices = meta_vars)
    updateSelectInput(session, "diff_covariates", choices = meta_vars)
    updateSelectInput(session, "pathway_strata", choices = c("None", meta_vars))
    updateSelectInput(session, "ml_outcome", choices = meta_vars)
    updateSelectInput(session, "ml_covariates", choices = meta_vars)
    updateSelectInput(session, "network_strata", choices = meta_vars)
    updateSelectInput(session, "biomarker_strata", choices = c("None", meta_vars))
    updateSelectInput(session, "meta_var_x", choices = meta_vars)
    updateSelectInput(session, "meta_var_y", choices = meta_vars)
    updateSelectInput(session, "meta_var_color", choices = meta_vars)
  }
}