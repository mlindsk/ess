#' Fit a decomposable graphical model
#' @description A generic method for structure learning in decomposable
#' graphical models
#' @param df Character data.frame
#' @param type Character ("fwd", "bwd", "tree" or "tfwd")
#' @param q Penalty term in the stopping criterion
#' where \code{0} = AIC and \code{1} = BIC. Anything in between is
#' referred to as \code{qic}
#' @param sparse_qic Logical. If \code{nrow(df)} is small, the tables
#' tends to be sparse. In these cases the usual penalty term of AIC and
#' BIC is often too restrictive. If \code{sparse_qic} is \code{TRUE}
#' this penality is computed according to a sparse criteria. The criteria
#' resembles the usual penalty as \code{nrow(df)} grows.
#' @param trace Logical indicating whether or not to trace the procedure
#' @param thres A threshold mechanism for choosing between two different ways of
#' calculating the entropy.
#' @param wrap logical specifying if the result of a run with type = "tree"
#' should be converted to a "fwd" object
#' @return A \code{gengraph} object representing a decomposable graph.
#' @examples
#'
#' g <- fit_graph(derma)
#' print(g)
#' plot(g)
#'
#' # Adjacency matrix and adjacency list
#' adjm <- adj_mat(g)
#' adjl <- adj_lst(g)
#' 
#' @details
#' The types are
#' \itemize{
#' \item "fwd": forward selection
#' \item "bwd": backward selection
#' \item "tree": Chow-Liu tree (first order interactions only)
#' \item "tfwd": A combination of "tree" and "fwd". This can speed up runtime considerably in high dimensions.
#' }
#' Using \code{adj_lst} on an object returned by \code{fit_graph} gives the
#' adjacency list corresponding to the graph. Similarly one can use \code{adj_mat}
#' to obtain an adjacency matrix. Applying the \code{rip} function on an
#' adjacency list returns the cliques and separators of the graph.
#' @references \url{https://arxiv.org/abs/1301.2267}, \doi{10.1109/ictai.2004.100} 
#' @seealso \code{\link{adj_lst}}, \code{\link{adj_mat}},
#' \code{\link{as_igraph}}, \code{\link{gengraph}}
#' @export
fit_graph <- function(df,
                      type       = "fwd",
                      q          = .5,
                      trace      = FALSE,
                      sparse_qic = FALSE,
                      thres      = 5,
                      wrap       = TRUE)
{

  n <- ncol(df)
  complete <- n * (n-1L) / 2L

  if (!is.data.frame(df)) stop("df must be a data.frame.")
  if (!(type %in% .types())) stop(.types_msg())
  if (q < 0) stop("q must be positive")
  
  x <- gengraph(df, type, q, sparse_qic)

  if (inherits(x, "fwd")) {
    if (!neq_empt_chr(as.vector(x$e))) {
      # If no edges are added in fwd_init, x$e = character(0)
      if (trace) msg(0L, complete, 0L, "delta-qic")
      return(x)
    }
  }  
  
  if (inherits(x, "tree")) return(fit_tree(x, df, wrap))
    
  triv     <- trivial(x, complete)
  update_k <- update_iteration(x)
  k        <- sum(x$adj_matrix)/2

  x <- walk(x = x, df = df, q = q, thres = thres)
  k <- update_k(k)

  if (k == triv) return(x)

  stp      <- stop_condition(x)
  stop_val <- attr(x$e, "d_qic")
  if (stp(stop_val)) return(x)

  while (!stp(stop_val)) {
    if (trace) msg(k, complete, stop_val, "delta-qic")
    x <- walk(x = x, df = df, q = q, thres = thres)
    k <- update_k(k)
    if (k == triv) {
      if (trace) msg(k, complete, stop_val, "delta-qic")
      return(x)
    }
    stop_val <- attr(x$e, "d_qic")
  }
  if (trace) msg(k, complete, stop_val, "delta-qic")
  return(x)
}


#' Fit a decomposable graphical model on each component
#' @description Structure learning in decomposable graphical models on
#' several components
#' @inheritParams fit_graph
#' @param comp A list with character vectors. Each element in the list is a
#' component in the graph (using expert knowledge)
#' @return An adjacency list object
#' @seealso \code{\link{fit_graph}}, \code{\link{adj_lst.gengraph}},
#' \code{\link{adj_mat.gengraph}}, \code{\link{walk.fwd}},
#' \code{\link{walk.bwd}}, \code{\link{gengraph}}
#' @export
fit_components <- function(df,
                      comp,
                      type   = "fwd",
                      q      = 0.5,
                      trace  = FALSE,
                      thres  = 5,
                      wrap   = TRUE)
{
  adj <- lapply(unname(comp), function(x) {
    fit_graph(df[, x, drop = FALSE], type = type,  q = q, trace = trace, thres = thres, wrap = wrap)
  })
  return(unlist(lapply(adj, adj_lst), recursive = FALSE))
}
