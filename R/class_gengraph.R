new_gengraph <- function(df, adj) {
  if (!setequal(colnames(df), names(adj))) stop("column names of df does not correspond to adj")
  structure(list(
    adj_list    = adj,                                         # graph as adjacency list
    adj_matrix  = as_adj_mat(adj),                             # graph as adjacency matrix
    adj_list_cg = NULL,                                        # clique list
    lvls        = .map_int(df, function(x) length(unique(x))), # level vector (for stopping criteria)
    mem         = list(ent = new.env(), npc = new.env())       # entropies and number of positive cells
  ),
  class = c("gengraph", "list")
  )
}

new_bwd <- function(df, q, sparse_qic) {
  adj  <- make_complete_graph(colnames(df))
  g    <- new_gengraph(df, adj)
  g$adj_list_cg <- rip(adj, check = FALSE)$C
  g$e  <- NULL # The newly deleted edge
  g$sparse_qic <- sparse_qic
  structure(g, class = c("bwd", class(g)))
}

new_fwd <- function(df, q, sparse_qic) {
  adj    <- make_null_graph(colnames(df))
  g      <- new_gengraph(df, adj)
  g$adj_matrix_cg <- as_adj_mat(make_complete_graph(colnames(df))) # Can be more efficient!
  g$adj_list_cg   <- as.list(names(adj))
  g$msi  <- list(S = NULL, max = list(e = character(0), idx = numeric(0), ins = vector("numeric", 2L)))
  g$e    <- new_edge()
  g$sparse_qic <- sparse_qic
  g      <- fwd_init(g, df, q)
  structure(g, class = c("fwd", class(g)))
}

new_tree <- function(df, sparse_qic) {
  adj    <- make_null_graph(colnames(df))
  g      <- new_gengraph(df, adj)
  g$adj_list_cg  <- as_adj_lst(adj)
  g$WGT  <- tree_weights(g, df) # Weights to use in kruskal procedure
  g$sparse_qic <- sparse_qic
  structure(g, class = c("tree", class(g)))
}

new_tfwd <- function(df, sparse_qic) {
  g <- fit_tree(new_tree(df), df, wrap = TRUE)
  g$sparse_qic <- sparse_qic
  structure(g, class = setdiff(c("tfwd", class(g)), "tree"))
}

new_edge <- function(e = character(0), d_qic = 0, idx = integer(0), ins = vector("integer", 2L)) {
  # e     : edge to be deletede or added
  # d_aic : entropy difference in the two competing models
  # idx   : in fwd procedure this is the index in msi where e lives
  # ins   : in fwd procedure this is the indicies in adj_list_cg where a new clique must be inserted
  structure(e, d_qic = d_qic, idx = idx, ins = ins)
}

#' A generic and extendable structure for decomposable graphical models
#' @description A generic structure for decomposable graphical models
#' @inheritParams fit_graph
#' @return A \code{gengraph} object with child class \code{type} used for model selection.
#' @examples
#'
#' gengraph(derma, type = "fwd")
#' gengraph(derma, type = "bwd")
#' 
#' @seealso \code{\link{adj_lst.gengraph}}, \code{\link{adj_mat.gengraph}}, \code{\link{fit_graph}}, \code{\link{walk.fwd}}, \code{\link{walk.bwd}}
#' @export
gengraph <- function(df, type = "fwd", q = .5, sparse_qic = TRUE) {
  switch(type,
    "fwd"  = new_fwd(df, q, sparse_qic),
    "bwd"  = new_bwd(df, q, sparse_qic),
    "tree" = new_tree(df, sparse_qic),
    "tfwd" = new_tfwd(df, sparse_qic)
  ) 
}

.types <- function() return(c("fwd", "bwd", "tree", "tfwd"))
.types_msg <- function() {
  paste0("Types must be in one of ",
    paste0(paste0(.types()[-length(.types())], collapse = ", "), " or ", .types()[length(.types())])
  )  
}
