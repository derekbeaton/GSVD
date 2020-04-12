#' is_GSVD
#'
#' Tests if the \code{x} object is of class type "GSVD"
#' @details The three primary functions in the \code{GSVD} package produce an inherited (hierarchical) class structure where all of them are of type "GSVD". Those functions are \code{\link{geigen}}, \code{\link{gsvd}}, and \code{\link{gplssvd}}.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class GSVD, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_GSVD <- function(x){
  inherits(x, "GSVD")
}

#' is_GSVD_geigen
#'
#' Tests if the \code{x} object is of class type "geigen"
#' @details Only \code{\link{geigen}} produces this class type.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class geigen, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_GSVD_geigen <- function(x){
  inherits(x, "geigen")
}

#' is_GSVD_gsvd
#'
#' Tests if the \code{x} object is of class type "gsvd"
#' @details Only \code{\link{gsvd}} produces this class type.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class gsvd, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_GSVD_gsvd <- function(x){
  inherits(x, "gsvd")
}

#' is_GSVD_gplssvd
#'
#' Tests if the \code{x} object is of class type "gplssvd"
#' @details Only \code{\link{gplssvd}} produces this class type.
#'
#' @param x object to test
#' @return boolean. \code{TRUE} if the object is of class gplssvd, FALSE otherwise.
#'
#' @seealso \code{\link{inherits}}
#'
#' @export
is_GSVD_gplssvd <- function(x){
  inherits(x, "gplssvd")
}



### strongly considering dropping this print method or radically changing it.
# #' @export
# print.GSVD <- function(x, ...){
#
#   if(!is_GSVD(x)){
#     stop("print.GSVD: Not a recognized class")
#   }
#
#   cat("**GSVD package object of class type 'GSVD'.**\n\n")
#   cat(sprintf("Number of total components = %d.\n", length(x$tau)))
#   cat(sprintf("Number of retained components = %d.\n", length(x$d)))
#
#   cat("The 'GSVD' object contains:\n")
#   gsvd_obj_mat <- cbind(
#     rbind("$d_full","Full set of singular values"),
#     rbind("$l_full","Full set of eigen values"),
#     rbind("$d","Retained set of singular values (k)"),
#     rbind("$l","Retained set of eigen values (k)")
#   )
#
#   print( format(gsvd_obj_mat, justify = "right"), quote="FALSE")
#   cat("**See specific print.*() method for more details on geigen(), gsvd(), and gplssvd() classes**\n")
#
# }


### strongly considering dropping this print method or radically changing it.
# #' @export
# summary.GSVD <- function(object, ...){
#   ## this inheritance is super dumb
#   x <- object
#   if(!is_GSVD(x)){
#     stop("summary.GSVD: Not a recognized class")
#   }
#
#   cat("**GSVD package object of class type 'GSVD'.**\n\n")
#   cat(sprintf("Number of total components = %d.\n", length(x$tau)))
#   cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))
#
#   print_mat <- round(cbind(x$tau[1:length(x$d)], cumsum(x$tau[1:length(x$d)]),  x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
#   print_mat <- format(print_mat, justify = "right")
#
#   colnames(print_mat) <- c("Explained variance","Cumulative explained variance","Eigenvalues", "Singular values")
#   rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
#   print(print_mat, quote=FALSE)
#
# }



#' @export
print.geigen <- function(x, ...){

  if(!is_GSVD(x)){
    stop("print.geigen: Not a GSVD object")
  }
  if(!is_GSVD_geigen(x)){
    stop("print.geigen: Not a recognized class")
  }

  cat("**GSVD package object of class type 'geigen'.**\n\n")
  cat(sprintf("geigen() was performed on a marix with %d columns/rows\n", nrow(x$v)))
  cat(sprintf("Number of total components = %d.\n", length(x$d_full)))
  cat(sprintf("Number of retained components = %d.\n", length(x$d)))

  cat("\nThe 'geigen' object contains:")
  geigen_obj_mat <- rbind(
    cbind("$d_full","Full set of singular values"),
    cbind("$l_full","Full set of eigen values"),
    cbind("$d","Retained set of singular values (k)"),
    cbind("$l","Retained set of eigen values (k)"),
    cbind("$v","Eigen/singular vectors"),
    cbind("$q","Generalized eigen/singular vectors"),
    cbind("$fj","Component scores")
  )
  rownames(geigen_obj_mat) <- rep("",nrow(geigen_obj_mat))
  colnames(geigen_obj_mat) <- rep("",ncol(geigen_obj_mat))
  print( geigen_obj_mat, quote="FALSE")

}

#' @export
summary.geigen <- function(object, ...){
  ## this inheritance is super dumb
  x <- object
  if(!is_GSVD(x)){
    stop("summary.geigen: Not a GSVD object")
  }
  if(!is_GSVD_geigen(x)){
    stop("summary.geigen: Not a recognized class")
  }

  cat("**GSVD package object of class type 'geigen'. Use 'print()' or 'print.geigen() for more information.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$d_full)))
  cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))

  print_mat <- round(cbind(x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
  print_mat <- format(print_mat, justify = "right")

  colnames(print_mat) <- c("Eigenvalues", "Singular values")
  rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
  print(print_mat, quote=FALSE)

}


#' @export
print.gsvd <- function(x, ...){

  if(!is_GSVD(x)){
    stop("print.gsvd: Not a GSVD object")
  }
  if(!is_GSVD_gsvd(x)){
    stop("print.gsvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gsvd'.**\n\n")
  cat("gsvd() was performed on a matrix with", nrow(x$u),"rows and", nrow(x$v),"columns\n")
  cat(sprintf("Number of components = %d.\n", length(x$d_full)))
  cat(sprintf("Number of retained components = %d.\n", length(x$d)))


  cat("\nThe 'gsvd' object contains:")
  gsvd_obj_mat <- rbind(
    cbind("$d_full","Full set of singular values"),
    cbind("$l_full","Full set of eigen values"),
    cbind("$d","Retained set of singular values (k)"),
    cbind("$l","Retained set of eigen values (k)"),
    cbind("$u","Left singular vectors (for rows of DAT)"),
    cbind("$v","Right singular vectors (for columns of DAT)"),
    cbind("$p","Left generalized singular vectors (for rows of DAT)"),
    cbind("$q","Right generalized singular vectors (for columns of DAT)"),
    cbind("$fi","Left component scores (for rows of DAT)"),
    cbind("$fj","Right component scores (for columns of DAT)")
  )
  rownames(gsvd_obj_mat) <- rep("",nrow(gsvd_obj_mat))
  colnames(gsvd_obj_mat) <- rep("",ncol(gsvd_obj_mat))
  print( gsvd_obj_mat, quote="FALSE")

}

#' @export
summary.gsvd <- function(object, ...){
  ## this inheritance is super dumb
  x <- object
  if(!is_GSVD(x)){
    stop("summary.gsvd: Not a GSVD object")
  }
  if(!is_GSVD_gsvd(x)){
    stop("summary.gsvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gsvd'. Use 'print()' or 'print.gsvd() for more information.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$d_full)))
  cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))

  print_mat <- round(cbind(x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
  print_mat <- format(print_mat, justify = "right")

  colnames(print_mat) <- c("Eigenvalues", "Singular values")
  rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
  print(print_mat, quote=FALSE)

}

#' @export
print.gplssvd <- function(x, ...){

  if(!is_GSVD(x)){
    stop("print.gplssvd: Not a GSVD object")
  }
  if(!is_GSVD_gplssvd(x)){
    stop("print.gplssvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gplssvd'.**\n\n")
  cat("gplssvd() was performed on an X matrix with", nrow(x$lx),"rows and", nrow(x$u),"columns and a Y matrix with" , nrow(x$ly),"rows and", nrow(x$v),"columns\n")
  cat(sprintf("Number of total components = %d.\n", length(x$d_full)))
  cat(sprintf("Number of retained components = %d.\n", length(x$d)))

  cat("\nThe 'gplssvd' object contains:")
  gplssvd_obj_mat <- rbind(
    cbind("$d_full","Full set of singular values"),
    cbind("$l_full","Full set of eigen values"),
    cbind("$d","Retained set of singular values (k)"),
    cbind("$l","Retained set of eigen values (k)"),
    cbind("$u","Left singular vectors (for columns of X)"),
    cbind("$v","Right singular vectors (for columns of Y)"),
    cbind("$p","Left generalized singular vectors (for columns of X)"),
    cbind("$q","Right generalized singular vectors (for columns of Y)"),
    cbind("$fi","Left component scores (for columns of X)"),
    cbind("$fj","Right component scores (for columns of Y)"),
    cbind("$lx","Left (X) latent variable scores (for rows of X)"),
    cbind("$ly","Right (Y) latent variable scores (for rows of Y)")
  )
  rownames(gplssvd_obj_mat) <- rep("",nrow(gplssvd_obj_mat))
  colnames(gplssvd_obj_mat) <- rep("",ncol(gplssvd_obj_mat))
  print( gplssvd_obj_mat, quote="FALSE")

}

#' @export
summary.gplssvd <- function(object, ...){

  ## this inheritance is super dumb
  x <- object
  if(!is_GSVD(x)){
    stop("summary.gplssvd: Not a GSVD object")
  }
  if(!is_GSVD_gplssvd(x)){
    stop("summary.gplssvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gplssvd'. Use 'print()' or 'print.gplssvd() for more information.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$d_full)))
  cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))

  print_mat <- round(cbind(x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
  print_mat <- format(print_mat, justify = "right")

  colnames(print_mat) <- c("Eigenvalues", "Singular values")
  rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
  print(print_mat, quote=FALSE)

}
