


copulareg <- function(y, x, var_type_y, var_type_x, ...) {
  UseMethod("copulareg")
}


copulareg <- function(y, x, var_type_y, var_type_x, ...) {

  copulareg.default(y, x, var_type_y, var_type_x, ...)

}


copulareg.default <- function(y, x, var_type_y, var_type_x,
                        family_set = c("gaussian", "clayton", "gumbel"),
                        extra_x = NULL, extra_y=NULL) {
  ##' copulareg
  ##' @aliases copulareg
  ##' @description This function fits joint distributions with an R-vine
  ##' pair copula structure, that is constructed in a specific way so that
  ##' the conditional density and distribution of the variable y can be
  ##' computed explicitly.
  ##' @param y A vector of n observations of the (univariate) outcome variable y
  ##' @param x A (n x p) matrix of n observations of p covariates
  ##' @param var_type_y A character that has to be specified as "d" or "c"
  ##' to indicate whether y is discrete or continuous, respectively.
  ##' @param var_type_x A vector of p characters that have to take the value
  ##' "c" or "d" to indicate whether each margin of the covariates is discrete
  ##' or continuous.
  ##' @param family_set A vector of strings that specifies the set of
  ##' pair-copula families that the fitting algorithm chooses from. For an
  ##' overview of which values that can be specified, see the documentation for
  ##' \link[rvinecopulib]{bicop}.
  ##' @param extra_x Optional extra values of x to use for estimating the
  ##' margins of the covariates.
  ##' @param extra_y Optional extra values of y to use to estimate the margin
  ##' of y.

  if (!is.null(extra_x) & !is.null(extra_y)) {
    fit <- .fit_model_xy(y, x, var_type_y, var_type_x, family_set,
                         .compute_distrbs(rbind(x, extra_x), var_type_x),
                         .compute_distrbs(c(y, extra_y), var_type_y))
  } else if (!is.null(extra_x) & is.null(extra_y)) {
    fit <- .fit_model_xy(y, x, var_type_y, var_type_x, family_set,
                         .compute_distrbs(rbind(x, extra_x), var_type_x))
  } else if (is.null(extra_x) & !is.null(extra_y)) {
    fit <- .fit_model_xy(y, x, var_type_y, var_type_x, family_set,
                         .compute_distrbs(c(y, extra_y), var_type_y))
  } else {
    fit <- .fit_model_xy(y, x, var_type_y, var_type_x, family_set)
  }

  class(fit) <- "copulareg"

  fit
}

.fit_model_xy <- function(y, x, var_type_y, var_types_x,
                         family_set=c("gaussian", "clayton", "gumbel"),
                         distr_x=NULL, distr_y=NULL) {
  # Fit a pair copula model, by first fitting a model to the p-dimensional X,
  # and then augmenting this with the variable y.

  if (is.null(distr_x)) {
    distr_x <- .compute_distrbs(x, var_types_x)
  }

  if (is.null(distr_y)) {
    distr_y <- .compute_distrbs(y, var_type_y)
  }

  model_x <- .fit_model_x(x, distr_x, var_types_x, family_set = family_set)

  .fit_model_y(y, distr_y, var_type_y, family_set, x, distr_x,
               model_x)
}

.fit_model_x <- function(x, distr_x, var_types_x, family_set) {
  # This function fits a pair copula model with an rvine structure
  # to the data x.
  u_x <- distr_x$transform(x, var_types_x)
  rvinecopulib::vinecop(data = u_x, var_types = var_types_x,
                        family_set = family_set)
}

.get_valid_edges <- function(t, rvine_mat, flipped_mat = F) {

  # This function finds all valid values that can
  # be entered into m_{d + 1 -t, 1} for t from 2 to d-1,
  # given an R-vine matrix

  if (t == 1) {
    seq(ncol(rvine_mat) - 1)
  } else {

    # Todo: This will always be false, so I should perhaps
    # rewrite the code below.
    if (!flipped_mat) {
      rvine_mat <- rvine_mat[seq(nrow(rvine_mat), 1), ]
    }

    k <- nrow(rvine_mat) + 2 - t
    match_vec <- rvine_mat[k:nrow(rvine_mat), 1]

    valid_values <- c()

    for (l in 2:(nrow(rvine_mat) + 1 - t)) {
      set_l_k <- rvine_mat[k:nrow(rvine_mat), l]
      if (setequal(set_l_k, match_vec)) {
        valid_values <- c(valid_values, rvine_mat[l, l])
      }
      if (k == nrow(rvine_mat)) {
        tilde_set_l_k <- rvine_mat[l, l]
      } else {
        tilde_set_l_k <- c(rvine_mat[l, l],
                           rvine_mat[(k + 1):nrow(rvine_mat), l])
      }
      if (setequal(tilde_set_l_k, match_vec)) {
        valid_values <- c(valid_values, rvine_mat[k, l])
      }
    }

    valid_values
  }
}

.get_hfunc <- function(j, cond_set, conditionals) {

  u <- get(.make_hfunc_key(j, cond_set),
           envir = conditionals)
  if (is.null(dim(u))) {
    u
  } else {
    u[, 1]
  }
}

.fit_model_y <- function(y, distr_y, var_type_y, family_set, x, distr_x,
                         model_x) {

  # Extract number of covariates, p
  p <- model_x$structure$d

  # Augment the matrix with a new zero row and column so that its dimension is
  # increased by 1 in both directions. The new variable (no. 4 in this case)
  # must be positioned in the bottom left corner, so this zero is replaced by
  # the new variable.
  #
  # Example:
  #    Original model matrix:
  #                          2 2 2
  #                          3 3 0
  #                          1 0 0
  #    New model matrix :
  #                          0 2 2 2
  #                          0 3 3 0
  #                          0 1 0 0
  #                          4 0 0 0


  new_mat <- cbind(c(rep(0, p), p + 1),
                   rbind(rvinecopulib::as_rvine_matrix(model_x$structure),
                         rep(0, p)))

  # Compute transformed distributions
  u_x <- distr_x$transform(x, model_x$var_types)

  # make a vinecopula object to store the (x, y) model in,
  # initially a copy of the model of x.
  model_xy <- model_x

  # Compute all the conditionals from the x - model
  conditionals <- .compute_conditionals(u_x, model_x)

  # Insert the pseudo-observations for y
  assign(.make_hfunc_key(p + 1, c()), distr_y$transform(y, var_type_y),
         envir = conditionals)

  for (t in 1:p) {

    # Compute which potential edges satisfy the proximity condition
    valid <- .get_valid_edges(t, new_mat)

    # Find the conditioning set, empty for t=1
    cond_set <- switch(1 + (t > 1), c(), new_mat[1:(t - 1), 1])

    if (length(valid) == 1) {
      # Only one possible edge
      new_mat[t, 1] <- valid
    } else {

      # More than one possibility, have to compute conditional Kendall's
      # coefficients, and select the variable that has the largest absolute
      # conditional Kendall's tau with the response.
      u_x <- sapply(valid, function(j) .get_hfunc(j, cond_set, conditionals))
      u_y <- .get_hfunc(p + 1, cond_set, conditionals)
      new_mat[t, 1] <- valid[which.max(abs(
        stats::cor(u_y, u_x, method = "kendall")
        ))]
    }

    # Combine the two variables that the new pair copula will be fit to.
    u_i_t <- .weave_transformed(get(.make_hfunc_key(p + 1, cond_set),
                                    envir = conditionals),
                                get(.make_hfunc_key(new_mat[t, 1],
                                                    cond_set),
                                    envir = conditionals))

    # Fit the new pair-copula
    if (t > length(model_xy$pair_copulas)) {
      model_xy$pair_copulas <- append(
        model_xy$pair_copulas,
        list(list(rvinecopulib::bicop(
          data = u_i_t,
          var_types = c(var_type_y,
                        model_x$var_types[new_mat[t, 1]]),
          family_set = family_set,
          par_method = "mle",
          selcrit = "aic"))), after = t)

    } else {
      model_xy$pair_copulas[[t]] <- append(
        model_xy$pair_copulas[[t]],
        list(rvinecopulib::bicop(
          data = u_i_t,
          var_types = c(var_type_y,
                        model_x$var_types[new_mat[t, 1]]),
          family_set = family_set,
          par_method = "mle",
          selcrit = "aic")), after = 0)
    }

    # update the log-likelihood
    model_xy$loglik <- model_xy$loglik + model_xy$pair_copulas[[t]][[1]]$loglik

    # Compute the new conditionals
    assign(.make_hfunc_key(p + 1, new_mat[1:t, 1]),
           .bicop_2_hbicop(u_i_t, bicop_obj = model_xy$pair_copulas[[t]][[1]],
                           cond_var = 2, return_u_minus = var_type_y == "d"),
           envir = conditionals)

  }

  # Complete the model object by setting the structure, the variable types,
  # and the total number of parameters
  model_xy$structure <- rvinecopulib::as_rvine_structure(new_mat)
  model_xy$var_types <- c(model_xy$var_types, var_type_y)
  model_xy$npars <- sum(
    sapply(1:(model_xy$structure$d - 1),
           function(t) sum(
             sapply(1:(model_xy$structure$d - t),
                    function(j) model_xy$pair_copulas[[t]][[j]]$npars))))

  list(model = model_xy, conditionals = conditionals, distr_x = distr_x,
       distr_y = distr_y, y = y)

}

.compute_distrbs <- function(x, x_type) {

  #
  # This function takes a matrix x, and returns a list containing empirical
  # cdf functions for each column, as well as a function that can transform
  # a new matrix that has the same number of columns, by these empirical cdf
  # functions.
  #

  # Check, in case x is 1d
  if (is.null(ncol(x))) {
    x <- as.matrix(x)
  }

  ecdf_np1 <- function(x, x_type) {
    #
    # copy of stats::ecdf, but with a different
    # dividing constant: (1/n) is replaced by (1/(n+1)), if x is continuous
    #
    x <- sort(x)
    n <- length(x)
    weight <- switch(x_type, c = 1 / (n + 1), d = 1 / n)
    if (n < 1)
      stop("'x' must have 1 or more non-missing values")
    vals <- unique(x)
    rval <- stats::approxfun(vals, cumsum(tabulate(match(x, vals))) * weight,
                             method = "constant", yleft = 0, yright = 1, f = 0,
                             ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
  }

  margins <- sapply(seq(ncol(x)), function(j) ecdf_np1(x[, j], x_type[j]))

  transform <- function(x, var_types_x) {
    # Check, in case x is 1d
    if (is.null(ncol(x))) {
      x <- as.matrix(x)
    }
    distr_plus <- sapply(seq_len(ncol(x)), function(i) margins[[i]](x[, i]))
    if (any(var_types_x == "d")) {
      distr_min <- sapply(which(var_types_x == "d"),
                          function(i) margins[[i]](x[, i] - 1))
      return(cbind(distr_plus, distr_min))
    } else{
      return(distr_plus)
    }
  }

  list(margins = margins, transform = transform)

}



.bicop_2_hbicop <- function(u, bicop_obj, cond_var = 2, return_u_minus = F) {

  # Function that lets you compute h-functions, without having to
  # specify the u_1^- when computing C(u_1 | u_2), when u_1 is
  # discrete. This is because u_1^- is redundant in that case, but rvinecopulibs
  # hbicop() function demands that it is provided. In addition, the specifying
  # return_u_minus = T, will make the function output the n x 2 matrix
  # [C(u_2 | u_1), C(u_2^- | u_1)] if cond_var = 1, or
  # [C(u_1 | u_2), C(u_1^- | u_2)] if cond_var = 2.

  if (!bicop_obj$var_types[c(2, 1)[cond_var]] == "d") {
    # In this case (that the conditioned variable is continuous), the
    # h-function in rvinecopulib can be called directly, as it then behaves as
    # expected

    if (return_u_minus) {
      cbind(rvinecopulib::hbicop(u, cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types),
            rvinecopulib::hbicop(u, cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types))
    } else {
      rvinecopulib::hbicop(u, cond_var = cond_var, family = bicop_obj$family,
                           rotation = bicop_obj$rotation,
                           parameters = bicop_obj$parameters,
                           var_types = bicop_obj$var_types)
    }
  } else {

    # This is more complicated. There are four_cases.

    u_columns_1 <- switch(1 + 2 * (cond_var - 1) + 1 * (ncol(u) == 4),
                          c(1, 2, 1, 3),
                          c(1, 2, 3, 4),
                          c(1, 2, 3, 2),
                          c(1, 2, 3, 4))

    if (return_u_minus) {

      u_columns_2 <- switch(1 + 2 * (cond_var - 1) + 1 * (ncol(u) == 4),
                            c(1, 3, 1, 3),
                            c(1, 4, 3, 4),
                            c(3, 2, 3, 2),
                            c(3, 2, 3, 4))

      cbind(rvinecopulib::hbicop(u[, u_columns_1], cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types),
            rvinecopulib::hbicop(u[, u_columns_2], cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types))

    } else {
      rvinecopulib::hbicop(u[, u_columns_1], cond_var = cond_var,
                           family = bicop_obj$family,
                           rotation = bicop_obj$rotation,
                           parameters = bicop_obj$parameters,
                           var_types = bicop_obj$var_types)
    }
  }
}

.make_hfunc_key <- function(i, cond_set) {

  #
  # This function creates a string key to be used when organising the evaluated
  # h-functions of every level of an R-vine pair copula model in a hash table,
  # using R's built-in enviroments.
  #

  cond_set <- sort(cond_set)
  paste0(c("(", i, "| ",
           paste0(c(sapply(cond_set[-length(cond_set)],
                           function(cond_var) paste0(c(cond_var, ", "),
                                              collapse = "")),
                    cond_set[length(cond_set)]),
                  collapse = ""),
           ")"), collapse = "")
}

.weave_transformed <- function(u1, u2) {

  # This function weaves two transformed variables (so that the resulting
  # variable has the form [u+, u-], if any of the two margins are discrete,
  # u- only contains columns from discrete variables),
  # where u_j will be a 2xn matrix if the conditioned variable
  # in the j-th transformed variable is discrete, and a numeric vector
  # of length n if continuous.

  u <- cbind(u1, u2)

  if (is.null(ncol(u1)) || ncol(u1) == 1) {
    u
  } else if (is.null(ncol(u2)) || ncol(u2) == 1) {
    u[, c(1, 3, 2)]
  } else {
    u[, c(1, 3, 2, 4)]
  }
}

.compute_conditionals <- function(u, model) {

  # Helper function that computes all of the transformed variables at every
  # level of the pair-copula model

  pick_u_cols <- function(var_inds, u, var_types) {
    # Small function that just selects the relevant columns of U for a specific
    # variable pair, used when computing the transformed variables for the
    # first tree, to avoid that the code in the loop below becomes too messy.
    if (any(var_types[var_inds] == "d")) {
      u[, c(var_inds, length(var_types) + cumsum(var_types == "d")[var_inds])]
    } else {
      u[, var_inds]
    }
  }

  # Create the hash table that we will store the transformed variables in,
  # and extract the R-vine matrix into a separate object to make the code in
  # the loop below easier to read.
  transformed_variables <- new.env(hash = TRUE)
  vinemat <- rvinecopulib::as_rvine_matrix(model$structure)
  d <- model$structure$d

  # First, add all the original variables
  for (j in 1:d) {
    assign(.make_hfunc_key(j, c()),
           pick_u_cols(c(j), u, model$var),
           envir = transformed_variables)
  }

  for (t in 1:(d - 1)) {
    for (e in 1:(d - t)) {
      # Extract the U-pair, and make the keys where we will insert the
      # new transformed variables we compute.
      if (t > 1) {
        # From the second tree, there are conditional variables D, stored in
        # the variable cond_set, and the variables U are new located in the
        # hash-table, so we'll have to extract them from there.
        cond_set <- vinemat[1:(t - 1), e]
        u_i_t <- .weave_transformed(get(.make_hfunc_key(vinemat[d + 1 - e, e],
                                                        cond_set),
                                        envir = transformed_variables),
                                    get(.make_hfunc_key(vinemat[t, e],
                                                        cond_set),
                                        envir = transformed_variables))

        key1 <- .make_hfunc_key(vinemat[d + 1 - e, e],
                                c(cond_set, vinemat[t, e]))
        key2 <- .make_hfunc_key(vinemat[t, e],
                                c(cond_set, vinemat[d + 1 - e, e]))
      } else {
        # In the second tree, there are no conditional variables D,
        # and the variables U are the empirical margins of the original
        # variables.
        u_i_t <- pick_u_cols(vinemat[c(d + 1 - e, t), e], u, model$var)
        key1 <- .make_hfunc_key(vinemat[d + 1 - e, e], vinemat[t, e])
        key2 <- .make_hfunc_key(vinemat[t, e], vinemat[d + 1 - e, e])
      }

      # Compute and insert the new transformed variables
      assign(key1,
             .bicop_2_hbicop(u_i_t,
                             model$pair_copulas[[t]][[e]], cond_var = 2,
                             return_u_minus = model$var[vinemat[d + 1 - e,
                                                                e]] == "d"),
             envir = transformed_variables)

      assign(key2,
             .bicop_2_hbicop(u_i_t,
                              model$pair_copulas[[t]][[e]], cond_var = 1,
                              return_u_minus = model$var[vinemat[t, e]] == "d"),
             envir = transformed_variables)
    }
  }
  transformed_variables
}
