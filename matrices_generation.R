library(igraph)

#' Function to generate a community graph (SBM logic)
#' @param n number of elements in the graph
#' @param prob list of probabilities of belonging to each group
#' @param p_in probability to have an edge between nodes from different communities
#' @param p_in probability to have an edge between nodes from the same community
community_graph <- function(n, prob = c(1/2,1/4,1/4), p_in = 0.5, p_out = 0.1,
                            block_sizes = NULL) {
  pref_mat <- matrix(p_out, length(prob), length(prob))
  diag(pref_mat) <- p_in
  if(is.null(block_sizes)) block_sizes<- c(rmultinom(1, n, prob))
  graph_mat <- as_adjacency_matrix(sample_sbm(n,
                                              pref.matrix = pref_mat,
                                              block.sizes = block_sizes ))
  graph_mat
}

#' Function to generate a possible value for Omega from a graph structure 
#' @param p_in probability to have an edge between nodes from different communities
generate_omega <- function(n, omega_structure = "community",
                           prob = c(1/2,1/4,1/4), p_in = 0.5, p_out = 0.1,
                           block_sizes = NULL, v = 0.3, u = 0.1){
  cond <- FALSE
  while(!cond){
    if(omega_structure == "community"){
      G <- community_graph(n, prob, p_in, p_out, block_sizes) 
    }
    
    # Ensuring that the network is not empty for AUC to make sense
    if(max(G) == 0){
      off_diag_indices <- which(row(matrix(1:n, n, n)) != col(matrix(1:n, n, n)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      G[selected_index[["row"]], selected_index[["col"]]] <- 1
      G[selected_index[["col"]], selected_index[["row"]]] <- 1
    }
    omega_tilde <- G * v
    omega <- omega_tilde + diag(abs(min(eigen(omega_tilde)$values)) + u, n, n)
    # Ensuring that the network is not full for AUC to make sense
    if(min(omega) > 0){ # Ensuring that the network has 0s for AUC to make sense
      off_diag_indices <- which(row(matrix(1:n, n, n)) != col(matrix(1:n, n, n)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      omega[selected_index[["row"]], selected_index[["col"]]] <- 0
      omega[selected_index[["col"]], selected_index[["row"]]] <- 0
    }
    cond <- all(eigen(omega)$values > 0)
  }
  as.matrix(omega)
}

# n = 20
# block_sizes = c(10, 5, 5)
# for(k in 1:10){
#   Om <- generate_omega(n, p_in = 0.9, p_out = 0.5,
#                        block_sizes = block_sizes, v = 0.3, u = 0.1)
#   write.csv(Om, paste0("matrices_1/om_before_break_", k, ".csv"), row.names=FALSE)
# }
# 
# for(k in 1:10){
#   Om <- generate_omega(n, p_in = 0.9, p_out = 0,
#                        block_sizes = block_sizes, v = 0.3, u = 0.1)
#   write.csv(Om, paste0("matrices_1/om_after_break_", k, ".csv"), row.names=FALSE)
# }


n = 20
block_sizes = c(10, 5, 5)
for(k in 1:10){
  A <- as.matrix(community_graph(n, p_in = 0.9, p_out = 0.5, block_sizes)) 
  write.csv(A, paste0("adjacency_matrices_1/A_before_break_", k, ".csv"), row.names=FALSE)
}

for(k in 1:10){
  A <- as.matrix(community_graph(n, p_in = 0.9, p_out = 0, block_sizes)) 
  write.csv(A, paste0("adjacency_matrices_1/A_after_break_", k, ".csv"), row.names=FALSE)
}
