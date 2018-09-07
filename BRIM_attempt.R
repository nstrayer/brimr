library(tidyverse)
library(igraph)
library(patchwork)
# set.seed(42)

# Helper functions
tr <- function(mat) sum(diag(mat))

# takes an n-by-c matrix
# returns a vector of length n with which of the c coluns was largest
maxRowEntry <- function(M){
  # x = M[1, ]
  apply(M, 1, function(x) which(x == max(x))[1]) %>% unlist()
} 

# Takes a vector of length n of assignments of cluster 1 to c
# Returns a matrix of dimension n-by-c
vecToIdMat <- function(vec, num_clust){
  # browser()
  purrr::map_df(vec, function(index){
    row <- rep(0, num_clust)
    row[index] <- 1L
    t(row) %>% 
      as_tibble() %>% {
        this <- .
        colnames(this) <- paste0('c', 1:num_clust)
        this
      }
  }) %>% 
    as.matrix()
}

# Takes a number of samples: size
# returns a vector with random assignment of cluster integer.
makeClusterIds <- function(size, num_clust){
  rep(1:num_clust, times = (size/num_clust) + 1 ) %>% 
    head(size)
  # sample(1:num_clust, size = size, replace = T) 
}

# takes two arrays of cluster assignments and an iteration number
# Returns a dataframe of the cluster assignment to node id plus the iteration number
collectClusterHist <- function(g1_assignments, g2_assignments, i){
  tibble(
    cluster = c(g1_assignments,g2_assignments),
    iteration = i
  ) %>% 
    mutate(id = 1:n())
}


network_edges <- read_csv('network_edges.csv')
num_clust <- 4

BRIM <- function(network_edges, num_clust){
  # basic settings
  max_num_steps <- 50
  sd_threshold <- 0.001
  sd_window <- 3
  num_steps <- 15
  
  adjacency_mat <- network_edges %>% 
    graph.data.frame() %>% 
    get.adjacency() %>% 
    as.matrix()
  
  m <- nrow(network_edges)
  n1 <- unique(network_edges$from) %>% length()
  n2 <- unique(network_edges$to) %>% length()
  n <- n1 + n2
  
  # indices of adjacency matrix for both node types
  g1_indices <- 1:n1
  g2_indices <- n1 + (1:n2)
  
  # Upper right of the adjacency matrix, aka the only really important part. 
  A_tilde <- adjacency_mat[g1_indices, g2_indices]
  
  # Matrix of null-model connection probabilities matching A_tilde (only need upper right block due to bipartite undirected nature.)
  # We want the number of expected edges to be equal to the number of edges observed in our network, thus prob of 1/m. 
  P_tilde <- matrix(1/m, nrow = n1, ncol = n2)
  
  B_tilde <- A_tilde - P_tilde

  # Calculate measure of clustering 'success' from the two cluster assignment matrices.
  getQ <- function(Rmat, Tmat) tr(t(Rmat) %*% B_tilde %*% Tmat)/m
  
  # cluster assignments for first node type...
  g1_assignments <- makeClusterIds(n1, num_clust)
  Rmat <- g1_assignments %>% vecToIdMat(num_clust)
  
  # ... and second node type
  g2_assignments <- makeClusterIds(n2, num_clust)
  Tmat <- g2_assignments %>% vecToIdMat(num_clust)
  
  Q_hist <- getQ(Rmat, Tmat)
  clust_hist <- collectClusterHist(g1_assignments, g2_assignments, i = 0)
  
  for (i in 1:max_num_steps) {
    # fix type 2 assignments and optimize type 1 assignments
    g1_assignments_new <- maxRowEntry(B_tilde %*% Tmat) %>% as.integer()
    Rmat <- vecToIdMat(g1_assignments_new, num_clust)
    
    # with new type 1 assignments, optimize the type 2 assignments
    g2_assignments <- maxRowEntry(t(B_tilde) %*% Rmat) %>% as.integer()
    Tmat <- vecToIdMat(g2_assignments, num_clust)
    
    # record the new Q value
    Q_hist <- c(Q_hist, getQ(Rmat, Tmat))
    
    # bundle the history's together
    current_clusters <- collectClusterHist(g1_assignments, g2_assignments, i)
    
    clust_hist <- bind_rows(
      clust_hist, 
      current_clusters
    )
    
    # check convergence termination threshold
    converged <- sd(tail(Q_hist, sd_window)) < sd_threshold
    if (converged) {
      print(sprintf('Converged in %i steps.', i))
      break
    }
  }
  
  Q_value_plot <- tibble(Q = Q_hist, iteration = as.numeric(0:i)) %>% 
    ggplot(aes(x = iteration, y = Q)) +
    geom_line() +
    scale_x_continuous(breaks = 0:i)

  cluster_size_plot <- clust_hist %>% 
    group_by(cluster, iteration) %>% 
    summarise(members = n()) %>% 
    ungroup() %>% 
    mutate(cluster = as.character(cluster)) %>% 
    ggplot(aes(x = iteration, y = members, fill = cluster)) +
    geom_col() +
    scale_x_continuous(breaks = 0:i) +
    labs(x = '')
  
  list(
    final_clusters = current_clusters %>% select(-i),
    cluster_history = clust_hist,
    cluster_history_plot = cluster_size_plot / Q_value_plot,
    number_of_iterations = i, 
    best_Q = tail(Q_hist, 1)
  )
}

Q_vals <- BRIM(network_edges, 2)$best_Q
for (c in 3:13) {
  Q_vals <- c(Q_vals, BRIM(network_edges, 10)$best_Q)
}


plot(Q_vals)
