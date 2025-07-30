# pathway_visualization.R

pathwayVisualizationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    box(
      title = "Pathway Visualization", width = 12, status = "success",
      div(
        id = ns("pathway_viewer"),
        style = "position: relative; width: 100%; height: 500px; overflow: auto;",
        imageOutput(ns("pathway_image"), height = "auto", width = "100%",
                    click = ns("image_click"))
      ),
      uiOutput(ns("zoom_controls")),
      htmlOutput(ns("mapped_data_info")),
      # Add debug info output
      verbatimTextOutput(ns("debug_info"))
    )
  )
}

pathwayVisualizationServer <- function(id, pathway_results, selected_pathway_index) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for image state
    pathway_img <- reactiveVal(NULL)
    img_zoom <- reactiveVal(1)
    img_pan <- reactiveVal(c(x = 0, y = 0))
    current_pathway_data <- reactiveVal(NULL)
    debug_info <- reactiveVal("Debug information will appear here")
    
    # Debug info output
    output$debug_info <- renderText({
      debug_info()
    })
    
    # Zoom controls UI
    output$zoom_controls <- renderUI({
      req(pathway_img())
      absolutePanel(
        top = 10, right = 10, fixed = TRUE,
        style = "background-color: rgba(255,255,255,0.8); padding: 5px; border-radius: 5px;",
        actionButton(ns("zoom_in"), "+", class = "btn-xs"),
        actionButton(ns("zoom_out"), "-", class = "btn-xs"),
        actionButton(ns("zoom_reset"), "Reset", class = "btn-xs")
      )
    })
    
    # Mapped data info
    output$mapped_data_info <- renderUI({
      req(current_pathway_data())
      
      pathway_data <- current_pathway_data()
      gene_ids <- tryCatch({
        unlist(strsplit(pathway_data$geneID, "/"))
      }, error = function(e) {
        debug_info(paste("Error splitting geneID:", e$message))
        return(NULL)
      })
      
      if(is.null(gene_ids)) return(NULL)
      
      n_genes <- length(gene_ids)
      mapped_genes <- if (!is.null(pathway_data$gene_data)) {
        sum(pathway_data$gene_data$Mapped == "Yes", na.rm = TRUE)
      } else {
        0
      }
      
      cpd_text <- if (!is.null(pathway_data$cpd_data)) {
        paste0("<br>Mapped compounds: ", nrow(pathway_data$cpd_data))
      } else {
        ""
      }
      
      HTML(paste0(
        "<div style='margin-top:10px; padding:10px; background-color:#f5f5f5; border-radius:5px;'>",
        "<strong>Mapped Data Summary</strong><br>",
        "Genes in pathway: ", n_genes, "<br>",
        "Mapped genes: ", mapped_genes, cpd_text,
        "</div>"
      ))
    })
    
    # Handle zoom controls
    observeEvent(input$zoom_in, {
      img_zoom(img_zoom() * 1.2)
    })
    
    observeEvent(input$zoom_out, {
      new_zoom <- img_zoom() / 1.2
      if (new_zoom >= 0.5) img_zoom(new_zoom)
    })
    
    observeEvent(input$zoom_reset, {
      img_zoom(1)
      img_pan(c(x = 0, y = 0))
    })
    
    # Handle panning
    observeEvent(input$image_click, {
      if (!is.null(input$image_click$clickId) && input$image_click$clickId == ns("pathway_image")) {
        current_pan <- img_pan()
        img_pan(c(x = input$image_click$x - current_pan["x"], 
                  y = input$image_click$y - current_pan["y"]))
      }
    })
    
    # Enhanced pathway visualization function with detailed debugging
    generate_pathway_visualization <- function(pathway_id, gene_list, cpd_list, species) {
      tryCatch({
        # Validate inputs
        if(is.null(pathway_id)) stop("Pathway ID is null")
        if(is.null(species)) stop("Species is null")
        
        # Create temp directory for pathview output
        temp_dir <- tempfile()
        dir.create(temp_dir)
        on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
        
        debug_info(paste("Temp directory created at:", temp_dir))
        
        # Remove species prefix if present
        pathway_id_clean <- sub(paste0("^", species), "", pathway_id)
        debug_info(paste("Cleaned pathway ID:", pathway_id_clean))
        
        # Generate a unique suffix
        vis_suffix <- paste0("vis_", format(Sys.time(), "%Y%m%d%H%M%S"))
        
        # Prepare data for pathview
        gene_data_debug <- if(!is.null(gene_list)) {
          paste("Gene data available with", length(gene_list), "entries")
        } else {
          "No gene data provided"
        }
        
        cpd_data_debug <- if(!is.null(cpd_list)) {
          paste("Compound data available with", length(cpd_list), "entries")
        } else {
          "No compound data provided"
        }
        
        debug_info(paste("Data status:", gene_data_debug, "|", cpd_data_debug))
        
        # Run pathview
        debug_info("Attempting to generate pathway visualization...")
        
        pv_out <- pathview::pathview(
          gene.data = gene_list,
          cpd.data = cpd_list,
          pathway.id = pathway_id_clean,
          species = species,
          kegg.native = TRUE,
          gene.idtype = "entrez",
          cpd.idtype = "kegg",
          out.suffix = vis_suffix,
          kegg.dir = temp_dir,
          low = list(gene = "blue", cpd = "blue"),
          mid = list(gene = "gray", cpd = "gray"),
          high = list(gene = "red", cpd = "red"),
          na.col = "white"
        )
        
        # Check output
        img_file <- file.path(temp_dir, paste0(species, pathway_id_clean, ".", vis_suffix, ".png"))
        
        if (file.exists(img_file)) {
          debug_info(paste("Successfully generated visualization at:", img_file))
          return(img_file)
        } else {
          debug_info("Pathview completed but no image file was created")
          return(NULL)
        }
      }, error = function(e) {
        debug_info(paste("Error in generate_pathway_visualization:", e$message))
        showNotification(paste("Pathway visualization error:", e$message), type = "error")
        return(NULL)
      })
    }
    
    # Generate pathway visualization when a pathway is selected
    observeEvent(selected_pathway_index(), {
      req(pathway_results(), selected_pathway_index())
      
      tryCatch({
        debug_info("Starting pathway visualization generation...")
        
        selected_pathway <- pathway_results()$kegg_res$ID[selected_pathway_index()]
        pathway_name <- pathway_results()$kegg_res$Description[selected_pathway_index()]
        
        debug_info(paste("Selected pathway:", selected_pathway, "-", pathway_name))
        
        # Generate visualization
        img_path <- generate_pathway_visualization(
          pathway_id = selected_pathway,
          gene_list = pathway_results()$gene_list,
          cpd_list = pathway_results()$cpd_list,
          species = pathway_results()$species
        )
        
        if (!is.null(img_path)) {
          pathway_img(img_path)
          img_zoom(1)
          img_pan(c(x = 0, y = 0))
          debug_info("Visualization successfully loaded")
          
          # Prepare mapped data for display
          pathway_df <- tryCatch({
            as.data.frame(pathway_results()$kegg_res[selected_pathway_index(), ])
          }, error = function(e) {
            debug_info(paste("Error getting pathway data:", e$message))
            return(NULL)
          })
          
          if(is.null(pathway_df)) return()
          
          # Extract mapped gene data
          gene_ids_in_pathway <- tryCatch({
            unlist(strsplit(pathway_df$geneID, "/"))
          }, error = function(e) {
            debug_info(paste("Error splitting geneID:", e$message))
            return(NULL)
          })
          
          if(is.null(gene_ids_in_pathway)) return()
          
          gene_data <- data.frame(
            GeneID = gene_ids_in_pathway,
            Log2FC = NA,
            Mapped = "No"
          )
          
          if (!is.null(pathway_results()$gene_list)) {
            matched_genes <- gene_ids_in_pathway[gene_ids_in_pathway %in% names(pathway_results()$gene_list)]
            if (length(matched_genes) > 0) {
              gene_data$Log2FC[match(matched_genes, gene_data$GeneID)] <- 
                pathway_results()$gene_list[matched_genes]
              gene_data$Mapped[match(matched_genes, gene_data$GeneID)] <- "Yes"
            }
          }
          
          # Extract mapped compound data
          cpd_data <- NULL
          if (!is.null(pathway_results()$cpd_list)) {
            cpd_data <- data.frame(
              CompoundID = names(pathway_results()$cpd_list),
              Log2FC = pathway_results()$cpd_list
            )
          }
          
          # Store current pathway data
          current_pathway_data(list(
            pathway_id = selected_pathway,
            pathway_name = pathway_name,
            gene_data = gene_data,
            cpd_data = cpd_data,
            geneID = pathway_df$geneID
          ))
          
        } else {
          debug_info("Visualization generation returned NULL")
          showNotification("Could not generate pathway visualization. Check debug info for details.", 
                           type = "warning")
        }
      }, error = function(e) {
        debug_info(paste("Error in visualization observer:", e$message))
        showNotification(paste("Visualization error:", e$message), type = "error")
      })
    })
    
    # Display pathway visualization with zoom/pan
    output$pathway_image <- renderImage({
      req(pathway_img())
      
      tryCatch({
        # Get image dimensions
        img_info <- magick::image_info(magick::image_read(pathway_img()))
        width <- img_info$width * img_zoom()
        height <- img_info$height * img_zoom()
        
        # Apply panning
        pan_x <- img_pan()["x"]
        pan_y <- img_pan()["y"]
        
        list(
          src = pathway_img(),
          contentType = "image/png",
          width = width,
          height = height,
          style = paste0("transform: translate(", pan_x, "px, ", pan_y, "px);"),
          alt = "Pathway visualization"
        )
      }, error = function(e) {
        debug_info(paste("Error rendering image:", e$message))
        showNotification("Error displaying pathway image", type = "error")
        return(list(src = ""))  # Return empty image on error
      })
    }, deleteFile = FALSE)
  })
}