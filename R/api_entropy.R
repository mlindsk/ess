na_tab <- function(df, a) {
  ct <- table(df[, a, drop = FALSE])
  names(dimnames(ct)) <- a # Needed for the onedimensional separators
  ct
}

joint_entropy <- function(df, npc = new.env()) {
  ## if( class(df) == "character" ) stop( "From entropy function: df is not a data.frame!" )
  x  <- na_tab(df, colnames(df))
  vx <- as.vector(x)
  nz <- vx[which(vx != 0)]
  Nx <- sum(nz)
  npc$value <- length(nz)
  -sum(nz/Nx * log(nz/Nx))
}

joint_entropy2 <- function(df, npc = new.env()) {
  A  <- apply(df, 1, paste0, collapse = ":")
  x  <- table(A)
  npc$value <- length(x)
  Nx <- sum(x)
  -sum(x/Nx * log(x/Nx))
}

#' Joint Entropy
#' 
#' @description Calculates the joint entropy over discrete variables in \code{df}
#' 
#' @param df data.frame
#' @param thres A threshold mechanism for choosing between two different ways of
#' calculating the entropy. Can Speed up the procedure with the "correct" value.
#' @param npc An environment. If supplied, the number of positive cells
#' in the underlying pmf will be stored in the environment with the name
#' \code{value}
#' @return A number representing the entropy of the variables in \code{df}.
#' @examples
#' entropy(derma[1:100, 1:3])
#' 
#' @export
entropy <- function(df, thres = 5, npc = new.env()) {
  if (ncol(df) <= thres) return(joint_entropy(df, npc))
  else return(joint_entropy2(df, npc))
}

entropy_difference <- function(e, S, df, mem, thres = 5) {
  # FIX: For large problems, delete some entropies?
  #      - maybe as for mem[["ent"]]ory accessible and
  #      - clean mem[["ent"]] if it reaches 50% of this or so?
  v <- unlist(es_to_vs(e))

  npc_S_val <- 0L
  H_S <- 0L
  if (neq_empt_chr(S)) {
    S_ <- sort_(S)
    if (exists(S_, envir = mem[["ent"]], inherits = FALSE)) {
      H_S <- mem[["ent"]][[S_]]
      npc_S_val <- mem[["npc"]][[S_]]
    } else {
      npc_S <- new.env()
      H_S   <- entropy(df[S], thres, npc_S)
      mem[["ent"]][[S_]] <- H_S
      mem[["npc"]][[S_]] <- npc_S[["value"]]
      npc_S_val <- npc_S[["value"]]
    }
  }

  npc_S_x_val <- 0L
  H_S_x <- 0L        
  Sx <- sort_(c(S, v[1]))
  if (exists(Sx, envir = mem[["ent"]], inherits = FALSE)) {
    H_S_x <- mem[["ent"]][[Sx]]
    npc_S_x_val <- mem[["npc"]][[Sx]]
  } else {
    npc_S_x <- new.env()
    H_S_x  <- entropy(df[c(S, v[1])], thres, npc_S_x)
    mem[["ent"]][[Sx]] <- H_S_x
    mem[["npc"]][[Sx]] <- npc_S_x[["value"]]
    npc_S_x_val <- npc_S_x[["value"]]
  }

  npc_S_y_val <- 0L
  H_S_y <- 0L
  Sy <- sort_(c(S, v[2]))
  if (exists(Sy, envir = mem[["ent"]], inherits = FALSE)) {
    H_S_y <- mem[["ent"]][[Sy]]
    npc_S_y_val <- mem[["npc"]][[Sy]]
  } else {
    npc_S_y <- new.env()
    H_S_y  <- entropy(df[c(S, v[2])], thres, npc_S_y)
    mem[["ent"]][[Sy]] <- H_S_y
    mem[["npc"]][[Sy]] <- npc_S_y[["value"]]
    npc_S_y_val <- npc_S_y[["value"]]
  }

  npc_S_xy_val <- 0L
  H_S_xy <- 0L
  Sxy <- sort_(c(S, v))
  if (exists(Sxy, envir = mem[["ent"]], inherits = FALSE)) {
    H_S_xy <- mem[["ent"]][[Sxy]]
    npc_S_xy_val <- mem[["npc"]][[Sxy]]
  } else {
    npc_S_xy <- new.env()
    H_S_xy  <- entropy(df[c(S, v)], thres, npc_S_xy)
    mem[["ent"]][[Sxy]] <- H_S_xy
    mem[["npc"]][[Sxy]] <- npc_S_xy[["value"]]
    npc_S_xy_val <- npc_S_xy[["value"]]
  }

  H_S_x_S_y <- H_S_x + H_S_y
  # Test needed to avoid < 0 due to floating point errors
  ent <- ifelse(isTRUE(all.equal(H_S_x_S_y, H_S_xy)), 0L,  H_S_x_S_y - H_S_xy - H_S)
  npc <- npc_S_x_val + npc_S_y_val - npc_S_xy_val - npc_S_val
  return(list(ent = ent, npc = npc))
}

