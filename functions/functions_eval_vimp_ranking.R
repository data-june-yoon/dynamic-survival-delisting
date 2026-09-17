
evaluate_mindepth_simulations <- function(depth_list, true_labels = c(1, 1, 1, 1, 0, 0, 0, 0)) {
    
    # Convert input list to a matrix
    if (is.list(depth_list) && !is.data.frame(depth_list)) {
        depth_matrix <- do.call(rbind, depth_list)
    } else {
        depth_matrix <- as.matrix(depth_list)
    }
    
    n_true <- sum(true_labels)
    n_noise <- length(true_labels) - n_true
    logical_truth <- as.logical(true_labels)
    
    # Calculate EVERY metric per row (per simulation)
    row_results <- t(apply(depth_matrix, 1, function(depth) {
        
        # --- Ratios (X1 < X5, etc.) ---
        # Evaluates to 1 (TRUE) or 0 (FALSE) for this specific row
        x1_ratio <- as.numeric(depth[1] < depth[5])
        x2_ratio <- as.numeric(depth[2] < depth[6])
        x3_ratio <- as.numeric(depth[3] < depth[7])
        x4_ratio <- as.numeric(depth[4] < depth[8])
        
        # --- Top-K Accuracy 
        # (= the proportion of your actual true variables that successfully landed in the top K ranked positions)
        ordered_indices <- order(depth, decreasing = FALSE)
        sorted_truth <- logical_truth[ordered_indices]
        top_k_accuracy <- sum(sorted_truth[1:n_true]) / n_true
        
        
        # Return all calculations as a named vector for this row
        c(x1_ratio = x1_ratio,
        x2_ratio = x2_ratio,
        x3_ratio = x3_ratio,
        x4_ratio = x4_ratio,
        Mean_Top_K_Acc = top_k_accuracy
        )
    }))
    
    final_summary <- as.data.frame(t(colMeans(row_results)))
    
    return(final_summary)
}

