## class methods and helpers

## need print and summary
  ### for now, I don't want any as.*() because I don't want others defining GSVD objects.




## is.*() methods

#' @export
is.GSVD <- function(x){
  inherits(x, "GSVD")
}
#' @export
is.geigen <- function(x){
  inherits(x, "geigen")
}
#' @export
is.gsvd <- function(x){
  inherits(x, "gsvd")
}
#' @export
is.gplssvd <- function(x){
  inherits(x, "gplssvd")
}


## GSVD package methods
## with "GSVD" I know that there will always be D, L, TAU; I should not worry about V, Q, and FJ for here.

#' @export
print.GSVD <- function(x){

  if(!is.GSVD(x)){
    stop("print.GSVD: Not a recognized class")
  }

  cat("**GSVD package object of class type 'GSVD'.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$tau)))
  cat(sprintf("Explained variance per component:\n", length(x$tau)))
  print( round(x$tau, digits = 3) )
  cat(sprintf("Number of retained components = %d.\n", length(x$d)))
  cat(sprintf("Cumulative variance for retained components = %f.\n", round(sum(x$tau[1:length(x$d)]), digits = 3)) )
  cat(sprintf("Explained variance per retained components:\n"))
  print( round(x$tau[1:length(x$d)], digits = 3) )


  cat("The 'GSVD' object contains:\n")
  gsvd_obj_mat <- cbind(
    rbind("$d.orig","Full set of singular values"),
    rbind("$l.orig","Full set of eigen values"),
    rbind("$tau","Percent explained variance per component for all components"),
    rbind("$d","Retained set of singular values (k)"),
    rbind("$l","Retained set of eigen values (k)")
  )

  print( format(gsvd_obj_mat, justify = "right"), quote="FALSE")
  cat("**See specific print.*() method for more details on geigen(), gsvd(), and gplssvd() classes**\n")

}

#' @export
summary.GSVD <- function(x){

  if(!is.GSVD(x)){
    stop("summary.GSVD: Not a recognized class")
  }

  cat("**GSVD package object of class type 'GSVD'.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$tau)))
  cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))

  print_mat <- round(cbind(x$tau[1:length(x$d)], cumsum(x$tau[1:length(x$d)]),  x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
  print_mat <- format(print_mat, justify = "right")

  colnames(print_mat) <- c("Explained variance","Cumulative explained variance","Eigenvalues", "Singular values")
  rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
  print(print_mat, quote=FALSE)

}



#' @export
print.geigen <- function(x){

  if(!is.geigen(x)){
    stop("print.geigen: Not a recognized class")
  }

  cat("**GSVD package object of class type 'geigen'.**\n\n")
  cat(sprintf("geigen() was performed on a marix with %d columns/rows\n", nrow(x$v)))
  cat(sprintf("Number of total components = %d.\n", length(x$tau)))
  cat(sprintf("Explained variance per component:\n", length(x$tau)))
  print( round(x$tau, digits = 3) )
  cat(sprintf("Number of retained components = %d.\n", length(x$d)))
  cat(sprintf("Cumulative variance for retained components = %f.\n", round(sum(x$tau[1:length(x$d)]), digits = 3)) )
  cat(sprintf("Explained variance per retained components:\n"))
  print( round(x$tau[1:length(x$d)], digits = 3) )

  cat("\nThe 'geigen' object contains:")
  geigen_obj_mat <- rbind(
    cbind("$d.orig","Full set of singular values"),
    cbind("$l.orig","Full set of eigen values"),
    cbind("$tau","Percent explained variance per component for all components"),
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
summary.geigen <- function(x){

  if(!is.geigen(x)){
    stop("summary.geigen: Not a recognized class")
  }

  cat("**GSVD package object of class type 'geigen'.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$tau)))
  cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))

  print_mat <- round(cbind(x$tau[1:length(x$d)], cumsum(x$tau[1:length(x$d)]),  x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
  print_mat <- format(print_mat, justify = "right")

  colnames(print_mat) <- c("Explained variance","Cumulative explained variance","Eigenvalues", "Singular values")
  rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
  print(print_mat, quote=FALSE)

}


#' @export
print.gsvd <- function(x){


  if(!is.gsvd(x)){
    stop("print.gsvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gsvd'.**\n\n")
  cat("gsvd() was performed on a matrix with", nrow(x$u),"rows and", nrow(x$v),"columns\n")
  cat(sprintf("Number of components = %d.\n", length(x$tau)))
  cat(sprintf("Explained variance per component:\n", length(x$tau)))
  print( round(x$tau, digits = 3) )
  cat(sprintf("Number of retained components = %d.\n", length(x$d)))
  cat(sprintf("Cumulative variance for retained components = %f.\n", round(sum(x$tau[1:length(x$d)]), digits = 3)) )
  cat(sprintf("Explained variance per retained components:\n"))
  print( round(x$tau[1:length(x$d)], digits = 3) )


  cat("\nThe 'gsvd' object contains:")
  gsvd_obj_mat <- rbind(
    cbind("$d.orig","Full set of singular values"),
    cbind("$l.orig","Full set of eigen values"),
    cbind("$tau","Percent explained variance per component for all components"),
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
summary.gsvd <- function(x){

  if(!is.gsvd(x)){
    stop("summary.gsvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gsvd'.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$tau)))
  cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))

  print_mat <- round(cbind(x$tau[1:length(x$d)], cumsum(x$tau[1:length(x$d)]),  x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
  print_mat <- format(print_mat, justify = "right")

  colnames(print_mat) <- c("Explained variance","Cumulative explained variance","Eigenvalues", "Singular values")
  rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
  print(print_mat, quote=FALSE)

}

#' @export
print.gplssvd <- function(x){

  if(!is.gplssvd(x)){
    stop("print.gplssvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gplssvd'.**\n\n")
  cat("gplssvd() was performed on an X matrix with", nrow(x$lx),"rows and", nrow(x$u),"columns and a Y matrix with" , nrow(x$ly),"rows and", nrow(x$v),"columns\n")
  cat(sprintf("Number of total components = %d.\n", length(x$tau)))
  cat(sprintf("Explained variance per component:\n", length(x$tau)))
  print( round(x$tau, digits = 3) )
  cat(sprintf("Number of retained components = %d.\n", length(x$d)))
  cat(sprintf("Cumulative variance for retained components = %f.\n", round(sum(x$tau[1:length(x$d)]), digits = 3)) )
  cat(sprintf("Explained variance per retained components:\n"))
  print( round(x$tau[1:length(x$d)], digits = 3) )

  cat("\nThe 'gplssvd' object contains:")
  gplssvd_obj_mat <- rbind(
    cbind("$d.orig","Full set of singular values"),
    cbind("$l.orig","Full set of eigen values"),
    cbind("$tau","Percent explained variance per component for all components"),
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
summary.gplssvd <- function(x){


  if(!is.gplssvd(x)){
    stop("summary.gplssvd: Not a recognized class")
  }

  cat("**GSVD package object of class type 'gplssvd'.**\n\n")
  cat(sprintf("Number of total components = %d.\n", length(x$tau)))
  cat(sprintf("Number of retained components = %d.\n\n", length(x$d)))

  print_mat <- round(cbind(x$tau[1:length(x$d)], cumsum(x$tau[1:length(x$d)]),  x$l[1:length(x$d)], x$d[1:length(x$d)]), digits = 3)
  print_mat <- format(print_mat, justify = "right")

  colnames(print_mat) <- c("Explained variance","Cumulative explained variance","Eigenvalues", "Singular values")
  rownames(print_mat) <- paste0("Component: ", 1:length(x$d))
  print(print_mat, quote=FALSE)

}
