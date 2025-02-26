TrokaChat_folder_structure <- function(main_directory) {
  # Define the folder names
  folders <- c("csv1", "DEG Output", "DEG Cleaned", "Noise Model", "cell_cell comm final",
               "Tensor Rank Determination", "Decomposed Tensor", "2 Condition Tensor")
  
  subfolders <- c("Heatmaps", "Cleaned Factors")
  
  # Create the main folder
  main_folder <- file.path(main_directory, "TrokaChat")
  if (!dir.exists(main_folder)) {
    dir.create(main_folder)
  }
  
  # Create the subdirectories in the main folder
  for (folder in folders) {
    subdirectory <- file.path(main_folder, folder)
    if (!dir.exists(subdirectory)) {
      dir.create(subdirectory)
    }
  }
  
  # Create the sub-subdirectories in the TrokaChat Charts folder
  for (subfolder in subfolders) {
    sub_subdirectory <- file.path(main_folder, "Decomposed Tensor", subfolder)
    if (!dir.exists(sub_subdirectory)) {
      dir.create(sub_subdirectory)
    }
  }
}

perform_communication_analysis <- function(seurat_obj, overall_list, n_permutations, output_dir, clusters_col, sample_col) {
  library(dplyr)
  library(pbapply)
  library(Seurat)
  library(parallel)
  library(progress)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract unique samples and clusters
  samples <- as.character(unique(seurat_obj@meta.data[[sample_col]]))
  clusters <- sort(as.numeric(as.character(unique(seurat_obj@meta.data[[clusters_col]]))))
  
  # Initialize a list to store the resulting data frames
  sample_dfs <- list()
  
  # Iterate over each sample
  for (sample in samples) {
    
    # Initialize an empty data frame to store results for this sample
    sample_df <- data.frame(gene = character(), cluster = character(), pct = numeric(), stringsAsFactors = FALSE)
    
    sample_indices <- which(FetchData(seurat_obj, vars = sample_col) == sample)
    sample_subset <- seurat_obj[, sample_indices]
    
    # Iterate over each cluster
    for (cluster in clusters) {
      # Subset the sample_subset for the current cluster
      cluster_indices <- which(FetchData(sample_subset, vars = clusters_col) == cluster)
      
      # Check if any cells are found for the current cluster
      if (length(cluster_indices) == 0) {
        cat(sprintf("No cells found in sample '%s' for cluster '%s'. Skipping this cluster.\n", sample, cluster))
        next  # Skip this cluster if no cells are found
      }
      
      cluster_subset <- sample_subset[, cluster_indices]
      
      # Proceed with calculations for the current cluster...
      cluster_df <- overall_list[[as.character(cluster)]]
      for (gene in unique(cluster_df$gene)) {
        expressing_cells <- sum(cluster_subset@assays$RNA@counts[gene, ] > 0)
        total_cells <- ncol(cluster_subset)
        pct_expressing <- (expressing_cells / total_cells) * 100
        sample_df <- rbind(sample_df, data.frame(gene = gene, cluster = cluster, pct = pct_expressing, stringsAsFactors = FALSE))
      }
    }
    
    # Add the sample data frame to the list
    sample_dfs[[sample]] <- sample_df
  }
  
  
  # Combine all sample data frames into one data frame
  combined_df <- bind_rows(sample_dfs, .id = "sample")
  
  # Convert to a character vector for easier indexing
  shortidents <- unique(combined_df$sample)
  
  # Function to copy Seurat object
  copy_object <- function(obj) {
    return(obj)
  }
  
  set.seed(12345)  # Set a master seed for reproducibility
  seeds <- sample(1:1e6, n_permutations)
  
  all_conditions <- list()
  
  # Extract unique genes of interest
  genes_of_interest <- unique(combined_df$gene)
  
  # Filter the Seurat object to include only the genes of interest
  object_new <- subset(seurat_obj, features = genes_of_interest)
  object_new@meta.data$TrokaChat_clusters_new <- object_new@meta.data[[clusters_col]]
  
  # Start timing the loop
  start_time <- Sys.time()
  
  # Permutation test
  for (c in seq_along(shortidents)) {
    cat("Starting condition ", shortidents[c], "\n")
    
    # Identify cells for the current sample
    condition_indices <- which(object_new@meta.data$sample == shortidents[c])
    
    # Check if any cells are present for this sample
    if (length(condition_indices) == 0) {
      cat(sprintf("No cells found for sample '%s'. Skipping this condition.\n", shortidents[c]))
      next  # Skip to the next sample
    }
    
    condition_object <- subset(object_new, cells = rownames(object_new@meta.data)[condition_indices])
    
    # Continue with the rest of your analysis...
    filtered_combined_df <- combined_df %>%
      filter(sample == shortidents[c])
    
    # Initialize holder2 with unique reference cluster and gene columns
    holder2 <- filtered_combined_df %>%
      select(cluster, gene) %>%
      distinct() %>%
      mutate(cluster = as.character(cluster), gene = as.character(gene))
    
    # Initialize holder3 as a list to store each permutation's result
    holder3 <- vector("list", n_permutations)
    
    holder3 <- mclapply(1:n_permutations, function(a) {
      set.seed(seeds[a])  # Set a unique seed for each iteration
      cat(sprintf("\nStarting permutation iteration: %d\n", a))
      
      # Ensure independent copy of condition_object for each iteration
      object_new_copy <- copy_object(condition_object)
      
      sample_indices <- which(object_new_copy@meta.data$sample == shortidents[c])
      cluster_counts_before <- table(object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices])
      
      cat("Before permutation:\n")
      print(head(object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices]))
      
      # Permute the labels
      permuted_labels <- sample(object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices])
      object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices] <- permuted_labels
      
      cat("After permutation:\n")
      print(head(object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices]))
      
      object_new_copy@meta.data[['sample_clus_TrokaChat']] <- ifelse(object_new_copy@meta.data$sample == shortidents[c], 
                                                                     paste0(object_new_copy@meta.data$sample, "_", object_new_copy@meta.data$TrokaChat_clusters_new), 
                                                                     paste0(object_new_copy@meta.data$sample, "_", object_new_copy@meta.data[[clusters]]))
      object_new_copy <- SetIdent(object_new_copy, value = 'sample_clus_TrokaChat')
      
      cluster_counts_after <- table(object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices])
      clusters_changed <- setdiff(names(cluster_counts_before), names(cluster_counts_after[cluster_counts_before == cluster_counts_after]))
      
      if (length(clusters_changed) == 0) {
        cat("All cluster counts are the same\n")
      } else {
        cat("Cluster counts are different in clusters:", paste(clusters_changed, collapse = ", "), "\n")
      }
      
      # Initialize output_df with the structure of holder2 and an additional column for the permutation
      output_df <- holder2
      output_df[[paste0("perm_", a)]] <- NA
      
      # Calculate the gene expression percentages for the permuted clusters
      for (cluster in unique(object_new_copy@meta.data$TrokaChat_clusters_new)) {
        
        # Subset the sample subset for the current cluster
        cluster_subset <- subset(object_new_copy, subset = TrokaChat_clusters_new == cluster)
        
        # Get the corresponding data frame from filtered_combined_df
        filtered_df <- filtered_combined_df %>%
          filter(cluster == cluster)
        
        # List of genes in the filtered_df
        gene_list <- unique(filtered_df$gene)
        
        # Subset the counts matrix to include only the relevant genes and cells
        cluster_counts <- cluster_subset@assays$RNA@counts[gene_list, ]
        
        # Calculate the number of expressing cells for each gene
        expressing_cells <- rowSums(cluster_counts > 0)
        
        # Calculate the total number of cells in the cluster
        total_cells <- ncol(cluster_counts)
        
        # Calculate the percentage of cells expressing each gene
        pct_expressing <- (expressing_cells / total_cells) * 100
        
        # Create a data frame to store the results for this cluster
        result_df <- data.frame(gene = gene_list, cluster = cluster, percentage = pct_expressing)
        
        # Update the result in the output data frame using dplyr
        output_df <- output_df %>%
          left_join(result_df, by = c("gene", "cluster")) %>%
          mutate(!!paste0("perm_", a) := coalesce(percentage, get(paste0("perm_", a)))) %>%
          select(-percentage)
      }
      
      return(output_df)
    }, mc.cores = detectCores() - 1)
    
    # Combine the results of all permutations into one data frame for the current sample
    holder3_combined <- Reduce(function(x, y) {
      merge(x, y, by = c("gene", "cluster"), all = TRUE)
    }, holder3)
    
    # Store the final holder3 for the current sample
    all_conditions[[shortidents[c]]] <- holder3_combined
    
    cat("Finished condition ", shortidents[c], "\n")
  }
  
  # End timing the loop
  end_time <- Sys.time()
  total_time <- end_time - start_time
  
  cat("Total time taken for the loop:", total_time, "\n")
  
  # Function to calculate mean and standard deviation of "perm_X" columns for each row
  calculate_mean_sd <- function(df) {
    # Identify the columns that match "perm_X"
    perm_cols <- grep("^perm_", names(df), value = TRUE)
    
    # Calculate the mean for each row
    df$mean_perm <- rowMeans(df[, perm_cols], na.rm = TRUE)
    
    # Calculate the standard deviation for each row
    df$sd_perm <- apply(df[, perm_cols], 1, sd, na.rm = TRUE)
    
    return(df)
  }
  
  # Apply the function to each data frame in the list
  all_conditions <- lapply(all_conditions, calculate_mean_sd)
  
  calculate_z_p_values_one_sided <- function(real_df, perm_df) {
    # Merge the real data with the permutation results
    merged_df <- merge(real_df, perm_df, by = c("gene", "cluster"))
    
    # Calculate Z score
    merged_df$z_score <- (merged_df$pct - merged_df$mean_perm) / merged_df$sd_perm
    
    # Calculate one-sided p-value
    merged_df$p_value <- 1 - pnorm(merged_df$z_score)
    
    return(merged_df)
  }
  
  adjust_p_values <- function(df_list) {
    # Combine all p-values into a single vector for adjustment
    all_p_values <- unlist(lapply(df_list, function(df) df$p_value))
    
    # Adjust p-values using the Benjamini-Hochberg procedure
    adjusted_p_values <- p.adjust(all_p_values, method = "BH")
    
    # Split the adjusted p-values back into the original list structure
    start <- 1
    for (i in seq_along(df_list)) {
      n <- nrow(df_list[[i]])
      df_list[[i]]$adjusted_p_value <- adjusted_p_values[start:(start + n - 1)]
      start <- start + n
    }
    
    return(df_list)
  }
  
  # Ensure the names of the lists match
  # if (!all(names(all_conditions) == names(overall_list))) {
  #   stop("The names of all_conditions and overall_list do not match.")
  # }
  
  # Calculate Z scores and p-values for each condition
  z_p_values <- mapply(calculate_z_p_values_one_sided, sample_dfs, all_conditions, SIMPLIFY = FALSE)
  
  # Adjust p-values
  z_p_values_adjusted <- adjust_p_values(z_p_values)
  
  # Reorder columns and remove perm_ columns
  z_p_values_adjusted <- lapply(z_p_values_adjusted, function(df) {
    perm_cols <- grep("^perm_", names(df), value = TRUE)
    non_perm_cols <- setdiff(names(df), perm_cols)
    df <- df[, c(non_perm_cols)]
    return(df)
  })
  
  # Save each condition's results to a CSV file
  for (condition in names(z_p_values_adjusted)) {
    write.csv(z_p_values_adjusted[[condition]], file = file.path(output_dir, paste0(condition, "_results.csv")), row.names = FALSE)
  }
  
  # Display the structure of the result
  str(z_p_values_adjusted)
  
  return(z_p_values_adjusted)
}

TrokaChat.DEG <- function(object, samples, control_condition, shortidents, filepath, export_folder_1,
                          clusters, sample_clus, cluster_range, sample_species,
                          cluster_pct_thresh, DefaultAssay) {
  ## Selecting Correct Ligand-Receptor (LR) Database for Species based on sample_species ----
  # Then, it loads the respective database and creates a vector of signaling genes.
  # Input: Seurat object, sample_species ("human" or "mouse")
  # Output: all_signaling_genes_final
  Idents(object) <- factor(x = Idents(object), levels = mixedsort(levels(object)))
  Idents(object)
  
  # Construct the file path based on sample_species
  if (sample_species == "human") {
    file_name <- "human_ALL_PATHWAYS_UPDATED.xlsx"
  } else if (sample_species == "mouse") {
    file_name <- "All Signaling Pathways.xlsx"
  } else {
    # Handle the case when sample_species does not match either condition
    stop("Invalid sample_species value. Please provide 'human' or 'mouse'.")
  }
  
  # Construct the full file path using the package's internal directory
  file_path <- system.file("ligand_receptor_db", file_name, package = "TrokaChat")
  if (file_path == "") {
    stop("Ligand-receptor database file not found in the package.")
  }
  
  # Read the Excel file
  lr_database <- read.xlsx(file_path)
  
  # Extract values from the 'ligand' column
  ligand_values <- as.character(lr_database$ligand)
  
  # Extract values from the 'subunit_' columns
  subunit_columns <- grep("^subunit_", names(lr_database), value = TRUE)
  subunit_values <- unlist(lr_database[subunit_columns])
  subunit_values <- as.character(subunit_values)
  
  # Combine the values into a single vector
  signaling_genes <- c(ligand_values, subunit_values)
  all_signaling_genes_final <- signaling_genes
  
  
  ## Generating Counts Table based on clusters and shortidents. ----
  # Creates a table that shows the number of cells per cluster for each condition.
  # Input: Seurat object, clusters, shortidents
  # Output: n_cells_bycondition_percluster
  object <- SetIdent(object, value = clusters)
  Idents(object) <- factor(x = Idents(object), levels = mixedsort(levels(object)))
  #Idents(object)
  n_cells_bycondition_percluster<-table(object@active.ident, object@meta.data$sample)
  n_cells_bycondition_percluster <- n_cells_bycondition_percluster[ , shortidents]
  print("#cells")
  print(n_cells_bycondition_percluster)
  
  
  ## Create List of Cluster Combinations for DEG analysis ----
  # Creates a list of all possible cluster combinations for subsequent DEG analysis
  # Input: shortidents, cluster_range, n_cells_bycondition_percluster
  # Output: total_list
  group <- paste0(shortidents, "_vs_", shortidents[1])
  
  #list_of_lists <- vector("list", length(shortidents))
  
  # INSERT HERE
  
  # Step 1: Calculate cutoff for each condition
  # Adjust the method to calculate cutoffs for each condition
  cutoffs <- sapply(shortidents, function(condition) {
    cells_in_condition <- n_cells_bycondition_percluster[, condition]
    if(cluster_pct_thresh > 0) {
      return(sum(cells_in_condition) * (cluster_pct_thresh / 100))
    } else {
      return(3) # Use 3 as the cutoff when cluster_pct_thresh is 0
    }
  }, USE.NAMES = TRUE)
  
  # Adjust to properly handle cluster numbering starting from 0
  valid_clusters_per_condition <- lapply(shortidents, function(condition) {
    cell_counts <- n_cells_bycondition_percluster[, condition]
    # Find clusters meeting the criteria
    valid_clusters_indices <- which(cell_counts >= max(cutoffs[condition], 3))
    
    # Adjust cluster numbering to start from 0 and assign names
    valid_clusters_names <- if (length(valid_clusters_indices) > 0) {
      # Subtract 1 to adjust for 0-based numbering
      paste0(condition, "_", valid_clusters_indices - 1)
    } else {
      character(0) # Return an empty character vector if no valid clusters
    }
    
    return(valid_clusters_names)
  })
  
  print(valid_clusters_per_condition)
  
  # Step 3: Find clusters valid across all conditions
  # Strip prefixes to get just the numeric part of cluster identifiers
  cluster_numbers_per_condition <- lapply(valid_clusters_per_condition, function(clusters) {
    as.integer(sub(".*_", "", clusters)) # Extracts the numeric part
  })
  
  # Find common cluster numbers across all conditions
  common_cluster_numbers <- Reduce(intersect, cluster_numbers_per_condition)
  
  
  print("LOOK")
  print(common_cluster_numbers)
  
  
  # Step 4: Generate comparisons for each condition
  # Initialize total_list for storing comparison lists for each condition
  total_list <- list()
  
  # Generate comparisons for each condition
  for (condition in shortidents) {
    condition_list <- list()
    
    # Rebuild valid cluster identifiers for the current condition based on common clusters
    valid_clusters_for_condition <- paste0(condition, "_", common_cluster_numbers)
    
    for (cluster_name in valid_clusters_for_condition) {
      # Determine valid comparisons
      if (condition == control_condition) {
        # For the control condition, exclude the current cluster from its comparisons
        valid_comparisons <- valid_clusters_for_condition[valid_clusters_for_condition != cluster_name]
      } else {
        # For experimental conditions, comparisons are against the control condition
        # Limit to common clusters, ensuring the comparison is meaningful across conditions
        control_clusters_for_comparison <- paste0(control_condition, "_", common_cluster_numbers)
        # Exclude the cluster if it matches the numeric part of the experimental cluster
        cluster_num <- sub(".*_", "", cluster_name)  # Extract numeric part
        control_cluster_to_exclude <- paste0(control_condition, "_", cluster_num)
        valid_comparisons <- control_clusters_for_comparison[control_clusters_for_comparison != control_cluster_to_exclude]
      }
      
      # Assign the list of valid comparisons to the condition list
      condition_list[[cluster_name]] <- valid_comparisons
    }
    
    # Assign the condition list to the total list
    total_list[[condition]] <- condition_list
  }
  
  
  
  print("START")
  print(total_list)
  print("END")
  
  saveRDS(total_list, file = paste0(filepath,export_folder_1,"/total_list.rds"))
  
  
  
  ##Initial DEG (Differential Expression Gene) Analysis ----
  # Executes initial DEG analysis by running FindMarkers() for each cluster combination in total_list
  # Input: Seurat object, total_list, all_signaling_genes_final
  # Output: Objects named *_TrokaChat1, *_TrokaChat2, *_TrokaChat3, etc. in the global environment
  rm(list=ls(pattern="_TrokaChat1"))
  rm(list=ls(pattern="_TrokaChat2"))
  rm(list=ls(pattern="_TrokaChat3"))
  
  object <- SetIdent(object, value = sample_clus)
  DefaultAssay(object) <- DefaultAssay
  features = {}
  for (i in 1:length(total_list)){
    for(a in 1:length(total_list[[i]]))
      assign(print(paste0(group[i],paste0(gsub("\\D+", "", names(total_list[[i]][a]))),"_TrokaChat",i)), FindMarkers(object, logfc.threshold = 0, only.pos = FALSE, ident.1 = print(paste0(names(total_list[[i]][a]))), ident.2 = print(total_list[[i]][[a]]), features=intersect(rownames(object), all_signaling_genes_final)), envir = environment())
  }
  
  
  
  
  ## For each cluster, only those DEGs that have a p-value less than 0.05 and are present in at least one of the conditions tested are considered significant. These genes are then used in subsequent DEG analysis. ----
  # Input: DEG results from the previous step
  # Output: TC_overall_list
  # Create the datalist
  datalist <- list()
  
  for (i in seq_along(shortidents)) {
    pattern_matches <- grep(paste0("*._TrokaChat",i,"$"), ls(environment()), value = TRUE)
    pattern <- mixedsort(pattern_matches)
    pattern_list <- do.call("list", mget(pattern))
    
    # str(pattern_list)  # Uncomment if you want to print the structure of pattern_list
    
    pattern <- gsub(paste0("_TrokaChat",i), '', pattern)
    pattern <- gsub('\\D', '', pattern)
    
    for (a in seq_along(pattern_list)) {
      if (nrow(pattern_list[[a]]) > 0) {
        pattern_list[[a]]$cluster <- pattern[a]
      }
    }
    
    datalist[[i]] <- do.call(rbind, pattern_list)
  }
  
  TC_all <- do.call(rbind, datalist)
  
  print(head(TC_all))
  # Clean up the gene names and filter
  TC_all$gene <- gsub("^.*\\.","", rownames(TC_all))
  TC_all <- TC_all[TC_all$p_val_adj < 0.05,]
  TC_all <- TC_all[TC_all$gene %in% all_signaling_genes_final, ]
  TC_all <- TC_all[order(TC_all$cluster),]
  
  # Split, remove duplicates, and sort
  TC_all_split <- split(TC_all, f = TC_all$cluster)
  TC_all_split <- lapply(TC_all_split, function(x) x[!duplicated(x[c("gene")]), ])
  TC_all_split <- TC_all_split[mixedsort(names(TC_all_split))]
  
  # Print the structure
  str(TC_all_split)
  
  TC_overall_list =  TC_all_split
  #INSERT HERE
  
  print("NEW: CHECK")
  str(TC_overall_list)
  print("Length of total list")
  print(length(total_list))
  print("this is the length of tc overall list")
  print(length(TC_overall_list))
  ## Conduct DEG With filtered list of genes ----
  # Run DEG analysis again but now using only significant genes from TC_overall_list.
  # Input: Seurat object, total_list, TC_overall_list
  # Output: Objects named *_TrokaChat_a, *_TrokaChat_b, *_TrokaChat_c, etc. in the global environment
  rm(list=ls(pattern="_TrokaChat_"))
  
  for (i in 1:length(total_list)){
    for (a in 1:length(TC_overall_list)){
      assign(print(paste0(group[i],paste0(gsub("\\D+", "", names(total_list[[i]][a]))),"_TrokaChat_",letters[i])), FindMarkers(object, logfc.threshold = 0,features = TC_overall_list[[a]]$gene, only.pos = FALSE, ident.1 = print(paste0(names(total_list[[i]][a]))), ident.2 = print(total_list[[i]][[a]])), envir = environment())
    }
  }
  
  saveRDS(TC_overall_list, file = paste0(filepath,export_folder_1,"/TC_overall_list.rds"))
  ## Assembles DEG data and pulls genes for final DEG anaylsis to recover pct.2 ----
  # Aggregates the DEG analysis results and prepares for the final DEG analysis.
  # Input: Results from the second DEG analysis
  # Output: dataframe_names_final
  rm(list=ls(pattern="dataframe_names_"))
  
  # Preallocate a list to hold your dataframes
  dataframe_names_list <- vector("list", length(total_list))
  Pattern2_list <- vector("list", length(total_list))
  
  # First loop
  for (i in seq_along(total_list)) {
    dataframe_names_list[[i]] <- mget(mixedsort(ls(pattern = paste0("_TrokaChat_", letters[i]))))
  }
  
  # Second loop
  for (i in seq_along(total_list)) {
    Pattern2_list[[i]] <- gsub(paste0("_TrokaChat_", letters[i]), '', names(dataframe_names_list[[i]]))
    Pattern2_list[[i]] <- gsub('\\D', '', Pattern2_list[[i]])
    perm2 <- Pattern2_list[[i]]
    var <- dataframe_names_list[[i]]
    
    for (a in seq_along(perm2)) {
      var[[a]]$cluster <- perm2[[a]]
      var[[a]]$gene <- row.names(var[[a]])
    }
    dataframe_names_list[[i]] <- var
  }
  
  
  
  ## Run Final DEG ----
  # Perform the final DEG analysis using the filtered list of genes.
  # Input: Seurat object, total_list, dataframe_names_list
  # Output: dataframe_names_new_list
  rm(list=ls(pattern="TrokaChatnew"))
  rm(list=ls(pattern="_new"))
  
  # Create a list to hold new dataframes
  dataframe_names_new_list <- vector("list", length(total_list))
  
  # Outer loop
  for (i in seq_along(total_list)){
    
    # Initialize a list to temporarily hold the FindMarkers results for each 'a'
    FindMarkers_results <- vector("list", length(dataframe_names_list[[i]]))
    
    # Inner loop
    for (a in seq_along(dataframe_names_list[[i]])){
      holder <- dataframe_names_list[[i]]
      
      # Create the new variable name
      new_var_name <- paste0(group[i], paste0(gsub("\\D+", "", names(total_list[[i]][a]))), "_TrokaChatnew_", letters[i])
      print(new_var_name)
      print(names(total_list[[i]][a]))
      print(names(total_list[[i]][-a]))
      
      # Call FindMarkers and store the result in the list
      FindMarkers_results[[a]] <- FindMarkers(
        object,
        features = holder[[a]]$gene,
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = FALSE,
        ident.1 = names(total_list[[i]][a]),
        ident.2 = names(total_list[[i]][-a])
      )
    }
    # Collect all FindMarkers_results for this 'i' and store in dataframe_names_new_list
    dataframe_names_new_list[[i]] <- FindMarkers_results
  }
  
  
  ## Add the pct.2 values from this analysis back to DEG output ----
  # Replaces the pct.2 column in dataframe_names_list with the one from dataframe_names_new_list
  # Input: dataframe_names_new_list, dataframe_names_list
  # Output: dataframe_names_final
  # Loop over each 'a' in total_list
  for (a in seq_along(total_list)) {
    
    # Loop over each 'i' in dataframe_names_new_list[[a]]
    for (i in seq_along(dataframe_names_new_list[[a]])) {
      temp <- dataframe_names_new_list[[a]]
      temp2 <- dataframe_names_list[[a]]
      
      # Order rows by row.names
      temp[[i]] <- temp[[i]][order(row.names(temp[[i]])), ]
      temp2[[i]] <- temp2[[i]][order(row.names(temp2[[i]])), ]
      
      # Replace pct.2 column in temp2 with the one from temp
      temp2[[i]]$pct.2 <- temp[[i]]$pct.2
      
      # Update dataframe_names_list
      dataframe_names_list[[a]] <- temp2
    }
  }
  
  # Initialize a list to hold the final dataframes
  dataframe_names_final <- vector("list", length(total_list))
  str(dataframe_names_list)
  # Loop over each 'i' in total_list
  for (i in seq_along(total_list)) {
    # Concatenate dataframes together and store the result in dataframe_names_final
    dataframe_names_final[[i]] <- do.call(rbind, dataframe_names_list[[i]])
  }
  
  
  ## Export the csv outputs ----
  # Writes the final DEG results to CSV files, one for each group
  # Input: dataframe_names_final, filepath, export_folder_1, group
  # Output: CSV files
  for(i in seq_along(dataframe_names_final)){
    write.csv(dataframe_names_final[[i]],file=paste0(filepath,export_folder_1,"/",paste0(group[i]),".csv"))
  }
  
}

TrokaChat.DEG.nulldist <- function(object,
                                   shortidents,
                                   control_condition,
                                   filepath,
                                   export_folder,
                                   clusters,
                                   n.perms,
                                   assay = "SCT") {
  # Set the default assay based on the user's choice
  DefaultAssay(object) <- assay
  
  export_folder_1 <- export_folder
  group <- paste0(shortidents, "_vs_", shortidents[1])
  
  ## Create counts table ----
  object <- SetIdent(object, value = paste0(clusters))
  n_cells_bycondition_percluster <- FetchData(object, vars = c("ident", "sample")) %>%
    dplyr::count(ident, sample) %>%
    tidyr::spread(ident, n)
  write.csv(n_cells_bycondition_percluster,
            file = file.path(filepath, export_folder_1, "counts.csv"),
            row.names = FALSE)
  
  total_list <- readRDS(file = paste0(filepath, "csv1/total_list.rds"))
  TC_overall_list <- readRDS(file = paste0(filepath, "csv1/TC_overall_list.rds"))
  holder2 <- list()
  object_new <- object
  object_new@meta.data$TrokaChat_clusters_new <- object_new@meta.data[[clusters]]
  
  ## Helper functions
  
  rename_and_merge_object <- function(object_new, shortidents) {
    object_sub <- subset(x = object_new, subset = sample == shortidents[1])
    object_sub2 <- subset(x = object_new, subset = sample == shortidents[1])
    object_sub <- RenameCells(object_sub, "cell1")
    object_sub@meta.data$sample <- paste0(shortidents[1], "_NEW")
    merge(object_sub, object_sub2)
  }
  
  findMarkersAndCreateDF <- function(object, total_list, group, c, TC_overall_list) {
    dataframe_names_fornull <- vector("list", length(TC_overall_list))
    unique_ident1 <- unique(sapply(total_list[[c]], names))
    unique_ident2_lengths <- integer(length(unique_ident1))
    for (i in seq_along(TC_overall_list)) {
      ident1 <- names(total_list[[c]][i])
      ident2 <- total_list[[c]][[i]]
      args <- list(features = TC_overall_list[[i]]$gene,
                   logfc.threshold = 0,
                   min.pct = 0,
                   only.pos = FALSE,
                   ident.1 = ident1,
                   ident.2 = ident2)
      if (assay == "SCT") {
        args$recorrect_umi <- FALSE
      }
      dataframe_names_fornull[[i]] <- do.call(FindMarkers, c(list(object), args))
      dataframe_names_fornull[[i]]$gene <- row.names(dataframe_names_fornull[[i]])
      dataframe_names_fornull[[i]]$cluster <- as.numeric(sub(".*_", "", ident1))
      unique_ident2_lengths[which(unique_ident1 == ident1)] <- length(ident2)
    }
    dataframe_names_fornull
  }
  
  FindMarkersagain <- function(total_list, object, c, TC_overall_list, dataframe_names_fornull) {
    dataframe_names_fornull_new <- vector("list", length(TC_overall_list))
    for (i in seq_along(TC_overall_list)) {
      ident1 <- names(total_list[[c]][i])
      ident2 <- names(total_list[[c]][-i])
      args <- list(features = dataframe_names_fornull[[i]]$gene,
                   logfc.threshold = 0,
                   min.pct = 0,
                   only.pos = FALSE,
                   ident.1 = ident1,
                   ident.2 = ident2)
      if (assay == "SCT") {
        args$recorrect_umi <- FALSE
      }
      result <- do.call(FindMarkers, c(list(object), args))
      result$gene <- row.names(result)
      result$cluster <- as.numeric(sub(".*_", "", ident1))
      dataframe_names_fornull_new[[i]] <- result
    }
    dataframe_names_fornull_new
  }
  
  process_data <- function(a, dataframe_names_fornull, dataframe_names_fornull_new) {
    for (i in seq_along(dataframe_names_fornull_new)) {
      df_fornull <- tibble::rownames_to_column(dataframe_names_fornull[[i]], "row_name")
      df_fornull_new <- tibble::rownames_to_column(dataframe_names_fornull_new[[i]], "row_name")
      if ("gene" %in% colnames(df_fornull)) {
        df_fornull <- df_fornull[order(df_fornull$gene), ]
        df_fornull_new <- df_fornull_new[order(df_fornull_new$gene), ]
      } else {
        df_fornull <- df_fornull[order(df_fornull$row_name), ]
        df_fornull_new <- df_fornull_new[order(df_fornull_new$row_name), ]
      }
      df_fornull_updated <- df_fornull %>%
        dplyr::left_join(df_fornull_new %>% dplyr::select(row_name, pct.2),
                         by = "row_name", suffix = c("", "_new")) %>%
        dplyr::mutate(pct.2 = dplyr::coalesce(pct.2_new, pct.2)) %>%
        dplyr::select(-pct.2_new) %>%
        tibble::column_to_rownames("row_name")
      if (!"gene" %in% colnames(df_fornull_updated)) {
        df_fornull_updated <- tibble::rownames_to_column(df_fornull_updated, "gene")
      }
      dataframe_names_fornull[[i]] <- df_fornull_updated
    }
    dataframe_names_fornull_final <- do.call(rbind, dataframe_names_fornull)
    if (!"gene" %in% colnames(dataframe_names_fornull_final)) {
      dataframe_names_fornull_final <- tibble::rownames_to_column(dataframe_names_fornull_final, "gene")
    }
    dataframe_names_fornull_final
  }
  
  process_holder <- function(export_folder, group_element) {
    holder3 <- data.table::rbindlist(holder2, use.names = TRUE, fill = TRUE)[, lapply(.SD, mean), by = .(gene, cluster)]
    holder3 <- as.data.frame(holder3)
    col_order <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster")
    holder3 <- holder3[, col_order]
    write.csv(holder3,
              file = file.path(filepath, export_folder, paste0("nulldistavgs_", group_element, ".csv")),
              row.names = FALSE)
  }
  
  copy_object <- function(obj) {
    obj
  }
  
  save_hdf5 <- function(permutation_list, h5_filepath, group) {
    if (file.exists(h5_filepath)) file.remove(h5_filepath)
    h5createFile(h5_filepath)
    for (i in seq_along(permutation_list)) {
      group_name <- paste0(group, "_permutation_", i)
      h5createGroup(h5_filepath, group_name)
      h5write(permutation_list[[i]], h5_filepath, paste0(group_name, "/data"))
    }
  }
  
  ## Main processing
  set.seed(12345)
  seeds <- sample(1:1e6, n.perms)
  
  for (c in seq_along(total_list)) {
    if (c == 1) {
      object_temp <- object_new
      object_final <- rename_and_merge_object(object_temp, shortidents)
      activeAssay <- DefaultAssay(object_final)
      if (activeAssay == "SCT") {
        object_final <- PrepSCTFindMarkers(object = object_final)
      }
      
      holder2 <- parallel::mclapply(1:n.perms, function(a) {
        tryCatch({
          set.seed(seeds[a])
          object_final_copy <- copy_object(object_final)
          sample_indices <- which(object_final_copy@meta.data$sample == paste0(control_condition, "_NEW"))
          object_final_copy@meta.data$TrokaChat_clusters_new[sample_indices] <-
            sample(object_final_copy@meta.data$TrokaChat_clusters_new[sample_indices])
          object_final_copy@meta.data[["sample_clus_TrokaChat"]] <- ifelse(
            object_final_copy@meta.data$sample == control_condition,
            paste0(object_final_copy@meta.data$sample, "_", object_final_copy@meta.data[[clusters]]),
            paste0(object_final_copy@meta.data$sample, "_", object_final_copy@meta.data$TrokaChat_clusters_new)
          )
          object_final_copy <- SetIdent(object_final_copy, value = "sample_clus_TrokaChat")
          if (assay != "SCT") {
            DefaultAssay(object_final_copy) <- assay
          }
          total_list_special <- total_list
          labels <- sub("_", "_NEW_", names(total_list_special[[1]]))
          names(total_list_special[[1]]) <- labels
          dataframe_names_fornull <- findMarkersAndCreateDF(object_final_copy,
                                                            total_list_special, group, c,
                                                            TC_overall_list = TC_overall_list)
          dataframe_names_fornull_new <- FindMarkersagain(total_list = total_list,
                                                          object = object_final_copy, c,
                                                          TC_overall_list = TC_overall_list,
                                                          dataframe_names_fornull = dataframe_names_fornull)
          process_data(a, dataframe_names_fornull, dataframe_names_fornull_new)
        }, error = function(e) { NULL })
      }, mc.cores = parallel::detectCores() - 1)
      
      holder2 <- Filter(function(x) !any(sapply(x, function(y) any(is.na(y)))), holder2)
      process_holder(export_folder_1, group[c])
      saveRDS(holder2,
              file = file.path(filepath, export_folder, paste0("holder2_perm_", group[c], ".rds")))
      h5_filepath <- file.path(filepath, export_folder, paste0("nulldist_", group[c], ".h5"))
      save_hdf5(holder2, h5_filepath, group[c])
      
    } else {
      rm("dataframe_names_fornull", "holder2", "dataframe_names_fornull_new")
      holder2 <- list()
      
      holder2 <- parallel::mclapply(1:n.perms, function(a) {
        tryCatch({
          set.seed(seeds[a])
          object_new_copy <- copy_object(object_new)
          sample_indices <- which(object_new_copy@meta.data$sample == paste0(shortidents[c]))
          object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices] <-
            sample(object_new_copy@meta.data$TrokaChat_clusters_new[sample_indices])
          object_new_copy@meta.data[["sample_clus_TrokaChat"]] <- ifelse(
            object_new_copy@meta.data$sample == shortidents[c],
            paste0(object_new_copy@meta.data$sample, "_", object_new_copy@meta.data$TrokaChat_clusters_new),
            paste0(object_new_copy@meta.data$sample, "_", object_new_copy@meta.data[[clusters]])
          )
          object_new_copy <- SetIdent(object_new_copy, value = "sample_clus_TrokaChat")
          DefaultAssay(object_new_copy) <- assay
          dataframe_names_fornull <- findMarkersAndCreateDF(object_new_copy, total_list, group, c,
                                                            TC_overall_list = TC_overall_list)
          dataframe_names_fornull_new <- FindMarkersagain(total_list = total_list,
                                                          object = object_new_copy, c,
                                                          TC_overall_list = TC_overall_list,
                                                          dataframe_names_fornull = dataframe_names_fornull)
          process_data(a, dataframe_names_fornull, dataframe_names_fornull_new)
        }, error = function(e) { NULL })
      }, mc.cores = parallel::detectCores() - 1)
      
      holder2 <- Filter(function(x) !any(sapply(x, function(y) any(is.na(y)))), holder2)
      process_holder(export_folder_1, group[c])
      saveRDS(holder2,
              file = file.path(filepath, export_folder, paste0("holder2_perm_", group[c], ".rds")))
      h5_filepath <- file.path(filepath, export_folder, paste0("nulldist_", group[c], ".h5"))
      save_hdf5(holder2, h5_filepath, group[c])
    }
  }
}

create_cluster_mapping_from_seurat <- function(seurat_obj, cluster_num_col, cluster_name_col) {
  mapping <- seurat_obj@meta.data %>%
    dplyr::select(all_of(cluster_num_col), all_of(cluster_name_col)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!cluster_num_col := as.numeric(as.character(!!sym(cluster_num_col)))) %>%  # Convert factor to numeric
    dplyr::rename(
      combined_cluster_number = !!sym(cluster_num_col),
      Source_name = !!sym(cluster_name_col)
    ) %>%
    dplyr::mutate(Target_name = Source_name)
  
  mapping$cluster_name <- mapping$Source_name
  
  return(mapping)
}

prepare_data <- function(file_path, mapping) {
  
  columns_to_remove <- c("Pathway Label", "source_target", "source_ligand", "target_receptor", 
                         "Communication Score", "Uniqueness Score", "logTF Communication Score", 
                         "adjusted_p_value", "Condition")
  # Load the CSV file
  data <- read_csv(file_path)
  # Remove specified columns
  data <- data[, setdiff(names(data), columns_to_remove)]
  
  # Map cluster names
  mapping_source <- mapping %>% select(combined_cluster_number, Source_name)
  mapping_target <- mapping %>% select(combined_cluster_number, Target_name)
  result <- data %>%
    left_join(mapping_source, by = c("Source" = "combined_cluster_number")) %>%
    rename(source_cluster_name = Source_name) %>%
    left_join(mapping_target, by = c("Target" = "combined_cluster_number")) %>%
    rename(target_cluster_name = Target_name) %>%
    select(-Source, -Target) %>%
    rename(Source = source_cluster_name, Target = target_cluster_name)
  
  return(result)
}

process_and_create_tensor <- function(result, tensor_data_file, tensor_dim_file) {
  
  # Helper function to count zeros and non-zeros in the tensor
  count_zeros_and_non_zeros <- function(tensor) {
    # Get the underlying array data
    tensor_data <- tensor@data
    
    # Count zeros and non-zeros
    total_elements <- prod(dim(tensor_data))
    zero_count <- sum(tensor_data == 0)
    non_zero_count <- total_elements - zero_count
    
    # Calculate percentages
    zero_percentage <- (zero_count / total_elements) * 100
    non_zero_percentage <- (non_zero_count / total_elements) * 100
    
    # Return results as a list
    return(list(
      total_elements = total_elements,
      zero_count = zero_count,
      non_zero_count = non_zero_count,
      zero_percentage = zero_percentage,
      non_zero_percentage = non_zero_percentage
    ))
  }
  
  # Scale residuals and prepare the data
  result_scaled <- result %>%
    mutate(Scaled_Residuals = rescale(Residuals, to = c(0, 1)))
  
  # Split into experimental and control datasets
  experimental <- result_scaled %>%
    mutate(Scaled_Residuals = Scaled_Residuals) %>%
    arrange(desc(Scaled_Residuals))
  
  control <- result_scaled %>%
    mutate(Scaled_Residuals = 1 - Scaled_Residuals) %>%
    arrange(desc(Scaled_Residuals))
  
  # Add condition labels and merge
  experimental$Condition <- "Experimental"
  control$Condition <- "Control"
  merged_data <- bind_rows(experimental, control)
  
  # Define dimensions for the tensor
  sources <- unique(merged_data$Source)
  targets <- unique(merged_data$Target)
  ligand_receptors <- unique(merged_data$ligand_receptor)
  conditions <- c("Experimental", "Control")
  
  mode_names <- list(
    Source = unique(merged_data$Source),
    Target = unique(merged_data$Target),
    LigandReceptor = unique(merged_data$ligand_receptor),
    Tools = conditions
  )
  
  # Initialize the tensor array
  tensor_array <- array(0, dim = c(length(sources), length(targets), length(ligand_receptors), length(conditions)))
  
  # Populate the tensor array
  for (i in 1:nrow(merged_data)) {
    s <- which(sources == merged_data$Source[i])
    t <- which(targets == merged_data$Target[i])
    lr <- which(ligand_receptors == merged_data$ligand_receptor[i])
    condition_index <- which(conditions == merged_data$Condition[i])
    tensor_array[s, t, lr, condition_index] <- merged_data$Scaled_Residuals[i]
  }
  
  # Create the tensor object and assign dimension names
  tensor <- as.tensor(tensor_array)
  dimnames(tensor@data) <- list(Source = sources, Target = targets, Ligand_Receptor = ligand_receptors, Condition = conditions)
  
  # Export tensor data to CSV
  write.csv(as.data.frame.table(tensor@data), file = tensor_data_file, row.names = FALSE)
  
  # Export tensor dimensions to CSV
  write.csv(data.frame(dimensions = dim(tensor@data)), file = tensor_dim_file, row.names = FALSE)
  
  # Calculate and print zero and non-zero counts
  counts <- count_zeros_and_non_zeros(tensor)
  cat("Total elements in tensor:", counts$total_elements, "\n")
  cat("Number of zeros:", counts$zero_count, "\n")
  cat("Number of non-zeros:", counts$non_zero_count, "\n")
  cat("Percentage of zeros:", round(counts$zero_percentage, 2), "%\n")
  cat("Percentage of non-zeros:", round(counts$non_zero_percentage, 2), "%\n")
  
  # Return both the tensor and mode_names in a list
  return(list(tensor = tensor, mode_names = mode_names))
}

load_and_map_factors <- function(dir_path, num_modes = 4, factor_prefix = "factor_", mapping_prefix = "mapping_") {
  # Initialize lists to store factors and mappings
  factors <- list()
  mappings <- list()
  
  # Load factor and mapping CSV files for each mode
  for (i in 0:(num_modes - 1)) {
    factor_file <- paste0(dir_path, "/", factor_prefix, i, ".csv")
    mapping_file <- paste0(dir_path, "/", mapping_prefix, i, ".csv")
    factors[[i + 1]] <- as.matrix(read.csv(factor_file))
    mappings[[i + 1]] <- read.csv(mapping_file)
  }
  
  # Function to map indices back to original categories
  map_to_categories <- function(factor, mapping) {
    # Assuming mapping has columns "Index" and "Category_0"
    rownames(factor) <- mapping$Category_0[match(0:(nrow(factor) - 1), mapping$Index)]
    return(factor)
  }
  
  # Apply mapping to each factor
  factors_mapped <- lapply(seq_along(factors), function(i) map_to_categories(factors[[i]], mappings[[i]]))
  
  # Create a CP-like result list
  cp_like_result <- list(U = factors_mapped)
  
  # Return both the mapped factors and the CP-like result
  return(list(factors_mapped = factors_mapped, cp_like_result = cp_like_result))
}

generate_and_export_heatmaps <- function(factors_mapped, mode_names, mode_labels, cluster_mapping, output_dir, heatmap_sizes = NULL) {
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define default heatmap sizes if not provided
  if (is.null(heatmap_sizes)) {
    heatmap_sizes <- list(
      Source = c(8, 6),
      Target = c(8, 6),
      LigandReceptor = c(8, 6),
      Tools = c(8, 6)
    )
  }
  
  # Helper function to create heatmap for a single mode
  create_clustered_heatmap_with_gini <- function(data, mode, mode_label, mode_items, cluster_mapping) {
    # Convert to matrix
    mat <- as.matrix(data)
    
    # Rename columns to "Factor 1", "Factor 2", etc.
    colnames(mat) <- paste("Factor", 1:ncol(mat))
    
    # Apply cluster mapping if it's a Source or Target mode
    if (mode %in% c("Source", "Target")) {
      name_map <- setNames(cluster_mapping$cluster_name, cluster_mapping$combined_cluster_number)
      new_row_names <- name_map[as.character(mode_items)]
      rownames(mat) <- ifelse(is.na(new_row_names), mode_items, new_row_names)
    } else {
      rownames(mat) <- mode_items
    }
    
    # Z-score normalization for each factor (column) for the heatmap
    mat_scaled <- scale(mat)
    
    # Calculate Gini coefficients for rows and columns using original data
    gini_rows <- apply(mat, 1, function(x) {
      x_pos <- x[x > 0]  # Only positive values
      if (length(x_pos) > 0) DescTools::Gini(x_pos) else 0
    })
    gini_columns <- apply(mat, 2, function(x) {
      x_pos <- x[x > 0]  # Only positive values
      if (length(x_pos) > 0) DescTools::Gini(x_pos) else 0
    })
    
    # Create the bar annotations for Gini coefficients
    row_anno <- rowAnnotation(
      " " = anno_barplot(gini_rows, gp = gpar(fill = "blue"), axis_param = list(labels = NULL))  # No axis labels
    )
    col_anno <- HeatmapAnnotation(
      " " = anno_barplot(gini_columns, gp = gpar(fill = "blue"), axis_param = list(labels = NULL)),  # No axis labels
      which = "column"
    )
    
    # Define the color function
    max_abs <- max(abs(mat_scaled), na.rm = TRUE)
    col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("blue", "white", "red"))
    
    # Create the heatmap with Gini coefficient bar graphs
    heatmap_plot <- Heatmap(mat_scaled,
                            name = paste("Loading Z-score"), 
                            col = col_fun,
                            cluster_rows = TRUE,
                            cluster_columns = TRUE,
                            show_row_names = TRUE,
                            show_column_names = TRUE,
                            column_names_rot = 45,  # Rotate column names
                            left_annotation = row_anno,  # Add row Gini bar graph
                            top_annotation = col_anno  # Add column Gini bar graph
    )
    
    return(heatmap_plot)
  }
  
  # Generate and export heatmaps for each mode
  modes <- names(mode_names)
  for (mode in modes) {
    # Create the heatmap
    heatmap_plot <- create_clustered_heatmap_with_gini(
      data = as.data.frame(factors_mapped[[which(names(mode_names) == mode)]]),
      mode = mode,
      mode_label = mode_labels[mode],
      mode_items = mode_names[[mode]],
      cluster_mapping = cluster_mapping
    )
    
    # Get dimensions for this mode
    dims <- heatmap_sizes[[mode]]
    if (is.null(dims)) {
      dims <- c(8, 6)  # Default size if not specified for this mode
    }
    
    # Set PDF file path
    pdf_file <- paste0(output_dir, mode, "_heatmap.pdf")
    
    # Open PDF device with specified width and height
    pdf(file = pdf_file, width = dims[1], height = dims[2])
    
    # Draw the heatmap
    draw(heatmap_plot)
    
    # Close the PDF device
    dev.off()
  }
}

process_ligand_receptor_factors <- function(factors_mapped, mode_names, species, output_dir, significant_data_file = NULL) {
  # Standardize species input and validate
  species <- tolower(species)
  if (!species %in% c("mouse", "human")) {
    stop("Species must be either 'mouse' or 'human'.")
  }
  
  # Print species for debugging
  cat("Processing for species:", species, "\n")
  
  # --- Step 1: Create ranked LigandReceptor factors workbook ---
  cat("Creating ranked LigandReceptor factors workbook...\n")
  
  lr_factors <- factors_mapped[[3]]
  lr_names <- mode_names$LigandReceptor
  wb <- createWorkbook()
  
  create_ranked_df <- function(factor_vector, lr_names) {
    df <- data.frame(LigandReceptor = lr_names, Value = factor_vector) %>%
      arrange(desc(Value)) %>%
      mutate(Rank = row_number())
    return(df)
  }
  
  for (i in 1:ncol(lr_factors)) {
    factor_df <- create_ranked_df(lr_factors[, i], lr_names)
    addWorksheet(wb, sheetName = paste("Factor", i))
    writeData(wb, sheet = i, x = factor_df)
    addStyle(wb, sheet = i, style = createStyle(textDecoration = "BOLD"), rows = 1, cols = 1:3)
    setColWidths(wb, sheet = i, cols = 1, widths = 50)
  }
  
  rankings_file <- file.path(output_dir, "LigandReceptor_Rankings.xlsx")
  saveWorkbook(wb, rankings_file, overwrite = TRUE)
  
  # --- Step 2: Clean the dataframes and save as CSV ---
  cat("Cleaning dataframes and saving as CSV...\n")
  
  sheet_names <- excel_sheets(rankings_file)
  factors_list <- lapply(sheet_names, function(sheet) read_excel(rankings_file, sheet = sheet))
  names(factors_list) <- sheet_names
  
  clean_dataframe <- function(df) {
    df <- as_tibble(df) %>%
      filter(Value != 0) %>%
      select(-Rank)
    return(df)
  }
  
  factors_list <- map(factors_list, clean_dataframe)
  cleaned_dir <- file.path(output_dir, "Cleaned_Factors", "cleaned_LigandReceptor_factors")
  dir.create(cleaned_dir, showWarnings = FALSE, recursive = TRUE)
  map2(factors_list, sheet_names, ~write.csv(.x, file = file.path(cleaned_dir, paste0(.y, "_cleaned.csv")), row.names = FALSE))
  
  # --- Step 3: Process significant_data_file ---
  if (!is.null(significant_data_file)) {
    cat("Processing significant_data_file...\n")
    
    sheet_names <- excel_sheets(significant_data_file)
    wb <- createWorkbook()
    
    # Define format functions for mouse and human species
    format_mouse <- function(lr_name) {
      tokens <- strsplit(lr_name, "_")[[1]]
      
      formatted_tokens <- sapply(tokens, function(token) {
        # Handle hyphenated tokens
        if (grepl("-", token)) {
          parts <- strsplit(token, "-")[[1]]
          if (length(parts) > 1) {
            # Keep first part as is, capitalize first letter of parts after hyphen
            first_part <- parts[1]
            rest_parts <- sapply(parts[-1], function(p) {
              paste0(toupper(substr(p, 1, 1)), tolower(substr(p, 2, nchar(p))))
            })
            return(paste(c(first_part, rest_parts), collapse = "-"))
          }
        }
        
        # Standard gene name formatting for mouse - capitalize first letter only for all genes
        return(paste0(toupper(substr(token, 1, 1)), tolower(substr(token, 2, nchar(token)))))
      })
      
      return(paste(formatted_tokens, collapse = "_"))
    }
    
    format_human <- function(lr_name) {
      return(toupper(lr_name))
    }
    
    for (sheet in sheet_names) {
      df <- read_xlsx(significant_data_file, sheet = sheet)
      
      if ("LigandReceptor" %in% colnames(df)) {
        # Apply case formatting based on species
        if (species == "mouse") {
          df$LigandReceptor <- sapply(df$LigandReceptor, format_mouse)
        } else if (species == "human") {
          df$LigandReceptor <- sapply(df$LigandReceptor, format_human)
        }
      }
      
      addWorksheet(wb, sheet)
      writeData(wb, sheet, df)
    }
    
    corrected_file <- file.path(cleaned_dir, "significant_data_corrected.xlsx")
    saveWorkbook(wb, corrected_file, overwrite = TRUE)
  } else {
    cat("significant_data_file is NULL. Skipping processing.\n")
  }
}

generate_factor_networks <- function(excel_file_path, species, kg_dir) {

  # Standardize species input
  species <- tolower(species)
  if (!species %in% c("human", "mouse")) {
    stop("Species must be either 'human' or 'mouse'.")
  }
  
  # Define the path to the species-specific KG files in the package
  pkg_kg_dir <- system.file("kg_base_files", species, package = "TrokaChatML")
  if (pkg_kg_dir == "") {
    stop(paste("KG files for species", species, "not found in the package."))
  }
  
  # Create the user-specified kg_dir if it doesn’t exist
  dir.create(kg_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Copy the species-specific KG directory to kg_dir
  copied_kg_dir <- file.path(kg_dir, species)
  fs::dir_copy(pkg_kg_dir, copied_kg_dir, overwrite = TRUE)
  message("Copied KG files from ", pkg_kg_dir, " to: ", copied_kg_dir)
  
  # Get the names of the sheets in the Excel file
  sheet_names <- excel_sheets(excel_file_path)
  
  # Create a list to store the data frames from each sheet
  factors_list <- lapply(sheet_names, function(sheet) {
    read_excel(excel_file_path, sheet = sheet)
  })
  names(factors_list) <- sheet_names
  
  # Function to create monoplex and bipartite networks for each factor
  create_networks <- function(df, factor_name) {
    # Create a subdirectory for this factor within kg_dir
    factor_dir <- file.path(kg_dir, factor_name)
    dir.create(factor_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1. Create monoplex network (self-loops without weights)
    monoplex <- df %>%
      dplyr::select(LigandReceptor) %>%
      dplyr::rename("Source" = "LigandReceptor") %>%
      dplyr::mutate(Target = Source)
    
    # Save monoplex network without headers
    write_tsv(monoplex, file.path(factor_dir, "factor_factor.tsv"), col_names = FALSE)
    
    # 2. Create bipartite network (with weights)
    bipartite <- df %>%
      mutate(
        Receptors = str_extract(LigandReceptor, "(?<=_).*"),
        Receptors = str_split(Receptors, "_")
      ) %>%
      unnest(Receptors) %>%
      dplyr::select(LigandReceptor, Receptors, Value) %>%
      dplyr::rename("Source" = "LigandReceptor", "Target" = "Receptors", "Weight" = "Value")
    
    # Save bipartite network without headers
    write_tsv(bipartite, file.path(factor_dir, "factor_gene.tsv"), col_names = FALSE)
  }
  
  # Apply the network creation function to each factor
  iwalk(factors_list, ~create_networks(.x, .y))
  
  # Function to remove header from a single TSV file
  remove_header <- function(file_path) {
    # Read the file with headers
    data <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
    # Write back without headers
    write_tsv(data, file_path, col_names = FALSE)
    cat("Removed header from:", file_path, "\n")
  }
  
  # Function to process all TSV files in a directory and its subdirectories
  process_directory <- function(dir_path) {
    # Find all TSV files recursively
    tsv_files <- list.files(dir_path, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
    if (length(tsv_files) == 0) {
      cat("No TSV files found in:", dir_path, "\n")
    } else {
      # Process each file to remove its header
      walk(tsv_files, remove_header)
      cat("Processed", length(tsv_files), "TSV files in total.\n")
    }
  }
  
  # Remove headers from all TSV files in kg_dir
  process_directory(kg_dir)
  
  # Return invisible NULL to indicate the function’s side effect
  invisible(NULL)
}

update_chembl_to_drug_names <- function(root_dir, cores = 1) {
  # Load required libraries
  library(webchem)
  library(dplyr)
  library(readr)
  library(progress)
  library(parallel)
  library(doParallel)
  
  # Function to convert a ChEMBL ID to a drug name
  chembl_to_drug_name <- function(chembl_id) {
    result <- chembl_query(chembl_id, resource = "molecule")
    if (!is.null(result[[1]]$pref_name)) {
      return(result[[1]]$pref_name)
    } else {
      return(chembl_id)
    }
  }
  
  # Function to find all unique ChEMBL IDs across files
  find_unique_chembl_ids <- function(root_dir) {
    unique_ids <- character(0)
    subdirs <- list.dirs(root_dir, full.names = TRUE, recursive = TRUE)
    
    for (subdir in subdirs) {
      file_path <- file.path(subdir, "multiplex_drug_filtered_updated.tsv")
      if (file.exists(file_path)) {
        data <- read_tsv(file_path, show_col_types = FALSE)
        unique_ids <- union(unique_ids, unique(data$node))
      }
    }
    
    return(unique_ids)
  }
  
  # Function to process ChEMBL IDs in parallel
  process_chembl_ids <- function(chembl_ids, cores) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    # Export the chembl_to_drug_name function to the workers
    clusterExport(cl, "chembl_to_drug_name", envir = environment())
    
    pb <- progress_bar$new(
      format = "[:bar] :percent Processed: :current/:total ChEMBL IDs",
      total = length(chembl_ids),
      clear = FALSE,
      width = 60
    )
    
    results <- foreach(id = chembl_ids, .combine = c, .packages = "webchem") %dopar% {
      result <- chembl_to_drug_name(id)
      pb$tick()
      setNames(result, id)
    }
    
    stopCluster(cl)
    return(results)
  }
  
  # Function to update a single file with the ChEMBL ID mapping
  update_file <- function(file_path, chembl_mapping) {
    data <- read_tsv(file_path, show_col_types = FALSE)
    data$node <- sapply(data$node, function(id) chembl_mapping[id])
    write_tsv(data, file_path)
    cat("Updated:", file_path, "\n")
  }
  
  # Main process
  cat("Finding unique ChEMBL IDs...\n")
  unique_ids <- find_unique_chembl_ids(root_dir)
  cat("Found", length(unique_ids), "unique ChEMBL IDs.\n")
  
  cat("Processing ChEMBL IDs...\n")
  chembl_mapping <- process_chembl_ids(unique_ids, cores)
  
  cat("Updating files...\n")
  subdirs <- list.dirs(root_dir, full.names = TRUE, recursive = TRUE)
  for (subdir in subdirs) {
    file_path <- file.path(subdir, "multiplex_drug_filtered_updated.tsv")
    if (file.exists(file_path)) {
      update_file(file_path, chembl_mapping)
    }
  }
  
  cat("All files have been processed.\n")
}

update_mondo_to_disease_names <- function(root_dir, cross_ref_file = "mondo_cross_ref/Cross-reference.txt") {
  # Load required libraries
  library(dplyr)
  library(readr)
  library(stringr)
  
  # Function to create MONDO ID to disease name mapping
  create_mondo_mapping <- function(cross_reference_path) {
    # Read the cross-reference data
    cross_reference <- read_delim(cross_reference_path, delim = "\t", show_col_types = FALSE)
    
    # Create mapping
    mondo_mapping <- cross_reference %>%
      select(MONDO, label) %>%
      separate_rows(MONDO, sep = "\\|") %>%
      separate_rows(label, sep = "\\|") %>%
      filter(!is.na(MONDO) & MONDO != "") %>%
      distinct()
    
    # Convert to named vector for faster lookup
    setNames(mondo_mapping$label, mondo_mapping$MONDO)
  }
  
  # Function to update a single file using the MONDO ID mapping
  update_file <- function(file_path, mondo_mapping) {
    # Read the TSV file
    data <- read_tsv(file_path, show_col_types = FALSE)
    
    # Update node column
    data$node <- sapply(data$node, function(id) {
      if (str_detect(id, "^MONDO_")) {
        return(mondo_mapping[id] %||% id)  # Use %||% for NULL coalescing
      } else {
        return(id)
      }
    })
    
    # Write updated data back to the file
    write_tsv(data, file_path)
    cat("Updated:", file_path, "\n")
  }
  
  # Main process
  cat("Creating MONDO ID to disease name mapping...\n")
  
  # Locate the cross-reference file within the package
  cross_reference_path <- system.file(cross_ref_file, package = "TrokaChatML")
  if (cross_reference_path == "") {
    stop("Cross-reference file not found in the package at: ", cross_ref_file)
  }
  
  # Create the mapping
  mondo_mapping <- create_mondo_mapping(cross_reference_path)
  cat("Mapping created with", length(mondo_mapping), "entries.\n")
  
  # Process all files
  cat("Updating files...\n")
  subdirs <- list.dirs(root_dir, full.names = TRUE, recursive = TRUE)
  for (subdir in subdirs) {
    file_path <- file.path(subdir, "multiplex_disease_filtered.tsv")
    if (file.exists(file_path)) {
      update_file(file_path, mondo_mapping)
    }
  }
  
  cat("All files have been processed.\n")
  invisible(NULL)
}