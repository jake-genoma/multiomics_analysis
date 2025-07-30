# modules/dataupload_helpers.R

# Helper function to validate and load data
load_omics_data <- function(file_input, header, sep, type, id_type, species = "human") {
  req(file_input)
  
  tryCatch({
    # Read data
    df <- read.delim(file_input$datapath, 
                     header = header, 
                     sep = sep,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
    
    # Basic validation
    if (ncol(df) < 2) {
      stop("Data must have at least 2 columns (identifiers + samples)")
    }
    
    # Check for duplicate identifiers
    if (any(duplicated(df[,1]))) {
      showNotification(paste("Duplicate identifiers found in", type, "data"), 
                       type = "warning")
    }
    
    # Check numeric values
    numeric_cols <- sapply(df[,-1, drop = FALSE], is.numeric)
    if (!all(numeric_cols)) {
      showNotification(paste("Non-numeric values detected in", type, "data. Converting where possible."), 
                       type = "warning")
      df[,-1] <- lapply(df[,-1, drop = FALSE], function(x) as.numeric(as.character(x)))
    }
    
    # Store original IDs for conversion tracking
    attr(df, "original_ids") <- df[,1]
    attr(df, "id_type") <- id_type
    attr(df, "species") <- species
    
    return(df)
  }, error = function(e) {
    showNotification(paste("Error loading", type, "data:", e$message), type = "error")
    return(NULL)
  })
}

# Get appropriate organism database
get_org_db <- function(species) {
  switch(species,
         "human" = "org.Hs.eg.db",
         "mouse" = "org.Mm.eg.db",
         "rat" = "org.Rn.eg.db",
         "org.Hs.eg.db") # default to human
}

# Convert IDs to symbols
convert_to_symbols <- function(df, aggregate = FALSE) {
  req(df)
  
  species <- attr(df, "species")
  org_db <- get_org_db(species)
  id_type <- attr(df, "id_type")
  
  tryCatch({
    # Convert to gene symbols
    id_map <- clusterProfiler::bitr(attr(df, "original_ids"),
                                    fromType = id_type,
                                    toType = "SYMBOL",
                                    OrgDb = org_db)
    
    # Check mapping success
    if (nrow(id_map) == 0) {
      stop("No identifiers could be mapped to gene symbols")
    }
    
    # Check for 1:1 mapping
    mapping_stats <- table(id_map[[id_type]])
    multi_mapped <- sum(mapping_stats > 1)
    
    if (multi_mapped > 0) {
      showNotification(paste(multi_mapped, "identifiers have multiple gene symbol mappings"), 
                       type = "warning")
    }
    
    # Merge with data
    merged <- merge(df, id_map, 
                    by.x = colnames(df)[1], 
                    by.y = id_type)
    
    # If aggregation is selected, average multiple mappings
    if (aggregate && multi_mapped > 0) {
      merged <- aggregate(merged[,-c(1, ncol(merged))], 
                          by = list(merged$SYMBOL), 
                          FUN = mean, na.rm = TRUE)
      colnames(merged)[1] <- "GeneSymbol"
    } else {
      # Keep all mappings but warn about duplicates
      if (any(duplicated(merged$SYMBOL))) {
        showNotification("Duplicate gene symbols exist - consider enabling aggregation", 
                         type = "warning")
      }
      merged <- merged[,c("SYMBOL", setdiff(colnames(merged), c(colnames(df)[1], "SYMBOL")))]
      colnames(merged)[1] <- "GeneSymbol"
    }
    
    # Add conversion attributes
    attr(merged, "original_ids") <- attr(df, "original_ids")
    attr(merged, "id_type") <- "SYMBOL"
    attr(merged, "species") <- species
    
    return(list(data = merged, stats = list(
      original_ids = length(attr(df, "original_ids")),
      mapped = nrow(id_map),
      unique_symbols = length(unique(id_map$SYMBOL)),
      multi_mapped = multi_mapped
    )))
    
  }, error = function(e) {
    stop(paste("Symbol conversion failed:", e$message))
  })
}