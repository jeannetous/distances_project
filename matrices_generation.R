library(igraph)

#' Function to generate a community graph (SBM logic)
#' @param Q number of clusters
#' @param prob list of probabilities of intra-community edge
#' @param p_in probability to have an edge between nodes from different communities
#' @param p_in probability to have an edge between nodes from the same community
community_graph <- function(Q, prob = c(1/2,1/4,1/4), p_in = 0.5, p_out = 0.1) {
  pref_mat <- matrix(p_out, length(prob), length(prob))
  diag(pref_mat) <- p_in
  graph_mat <- as_adjacency_matrix(sample_sbm(Q,
                                              pref.matrix = pref_mat,
                                              block.sizes = c(rmultinom(1, Q, prob)) ))
  graph_mat
}

#' Function to generate a possible value for Omega from a graph structure 
#' @param p_in probability to have an edge between nodes from different communities
generate_omega <- function(Q, omega_structure = "community",
                           prob = c(1/2,1/4,1/4), v = 0.3, u = 0.1){
  cond <- FALSE
  while(!cond){
    # if(omega_structure == "erdos_reyni") G <- erdos_reyni_graph(Q)
    # if(omega_structure == "preferential_attachment") G <- preferential_attachment_graph(Q)
    if(omega_structure == "community") G <- community_graph(Q, prob)
    
    # Ensuring that the network is not empty for AUC to make sense
    if(max(G) == 0){
      off_diag_indices <- which(row(matrix(1:Q, Q, Q)) != col(matrix(1:Q, Q, Q)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      G[selected_index[["row"]], selected_index[["col"]]] <- 1
      G[selected_index[["col"]], selected_index[["row"]]] <- 1
    }
    omega_tilde <- G * v
    omega <- omega_tilde + diag(abs(min(eigen(omega_tilde)$values)) + u, Q, Q)
    # Ensuring that the network is not full for AUC to make sense
    if(min(omega) > 0){ # Ensuring that the network has 0s for AUC to make sense
      off_diag_indices <- which(row(matrix(1:Q, Q, Q)) != col(matrix(1:Q, Q, Q)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      omega[selected_index[["row"]], selected_index[["col"]]] <- 0
      omega[selected_index[["col"]], selected_index[["row"]]] <- 0
    }
    cond <- all(eigen(omega)$values > 0)
  }
  as.matrix(omega)
}


for(k in 1:10){
  Om <- generate_omega(10, prob = c(2/3,1/3,1/3))
  write.csv(Om, paste0("matrices_1/om_before_break_", k, ".csv"), row.names=FALSE)
}

for(k in 1:10){
  Om <- generate_omega(10, prob = c(1/2, 0, 0))
  write.csv(Om, paste0("matrices_1/om_after_break_", k, ".csv"), row.names=FALSE)
}
