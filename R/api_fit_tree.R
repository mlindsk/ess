tree_weights <- function(x, df) {
  # x: gengraph object
  nodes   <- colnames(df) 
  n       <- length(nodes) 
  pairs   <- utils::combn(nodes, 2,  simplify = FALSE) 
  weights <- structure(vector(mode = "numeric", length = n * (n - 1) / 2), names = "")
  for (j in 1:n) {
    npc_j <- new.env()
    x$mem[["ent"]][[nodes[j]]] <- entropy(df[nodes[j]], npc = npc_j)
    x$mem[["npc"]][[nodes[j]]] <- npc_j[["value"]]
  }
  for (p in seq_along(pairs)) {
    edge <- sort_(pairs[[p]])
    ed   <- entropy_difference(edge, character(0), df, x$mem)$ent
    weights[p] <- ed
    names(weights)[p] <- edge
  }
  
  return(sort(weights, decreasing = TRUE))
}


kruskal      <- function(x) UseMethod("kruskal")

kruskal.tree <- function(x) {
  n          <- length(x$adj_list)
  nodes      <- names(x$adj_list)
  node_pairs <- es_to_vs(names(x$WGT))
  number_of_nodes_total <- n
  number_of_nodes_added <- 0L
  for (e in seq_along(x$WGT)) {
    if (number_of_nodes_added == number_of_nodes_total - 1) return(x)
    node1 <- node_pairs[[e]][1]
    node2 <- node_pairs[[e]][2]
    component1 <- dfs(x$adj_list, node1)
    component2 <- dfs(x$adj_list, node2)
    if (!neq_empt_chr(intersect(component1, component2))) {
      x$adj_list[[node1]] <- c(x$adj_list[[node1]], node2)
      x$adj_list[[node2]] <- c(x$adj_list[[node2]], node1)
      x$adj_matrix[node1, node2] <- 1L                     
      x$adj_matrix[node2, node1] <- 1L                     
      number_of_nodes_added <- number_of_nodes_added + 1L
    }
  }
  # FIX: Update adj_list_cg - the print method is not correct
  return(x)
}


## as_fwd <- function(adj, mat, lst, ...) UseMethod("as_fwd")

tree_as_fwd <- function(x, df) {
  x$adj_list_cg <- rip(x$adj_list, check = FALSE)$C
  x$e    <- new_edge()
  nC     <- length(x$adj_list_cg)
  x$adj_matrix_cg <- matrix(0L, nC, nC)
  msi    <- vector("list", 0L)
  k         <- 1L
  if (nC > 1) {
    for (i in 2:nC) {
      for (j in 1:(i-1)) {
        Ci   <- x$adj_list_cg[[i]]
        Cj   <- x$adj_list_cg[[j]]
        Sij  <- intersect(Ci, Cj)
        if (neq_empt_chr(Sij)) { # Note: This ONLY work for trees
          x$adj_matrix_cg[i,j]  <- 1L
          x$adj_matrix_cg[j,i]  <- 1L
          Ci_minus_Sij <- setdiff(Ci, Sij)
          Cj_minus_Sij <- setdiff(Cj, Sij)
          edge_ij      <- sort_(c(Ci_minus_Sij, Cj_minus_Sij))
          ent_ij       <- entropy_difference(edge_ij, Sij, df, x$mem)$ent
          if (ent_ij >= attr(x$e, "d_qic")) {
            x$e <- new_edge(edge_ij, ent_ij, k, c(i, j))
          }
          msi[[k]] <- list(C1 = Ci, C2 = Cj, S = Sij, e = structure(ent_ij, names = edge_ij))
          k <- k + 1L
        }
      }
    }
  }
  x$msi <- msi
  class(x) <- setdiff(c("fwd", class(x)), "tree")
  return(x)
}

fit_tree <- function(x, df, wrap = TRUE) {
  if (wrap) return(tree_as_fwd(kruskal(x), df))
  else return(kruskal(x))
}
