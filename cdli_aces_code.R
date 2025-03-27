# define required packages
required_packages <- c("remotes", "aces", "readxl", "dplyr", "lpSolveAPI")

# install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "aces") {
      remotes::install_github("Victor-Espana/aces")
    } else {
      install.packages(pkg)
    }
  }
}

# load packages
library("aces")
library("readxl")
library("dplyr")
library("lpSolveAPI")

# ============================= #
# Directional Distance Function #
# ============================= #

# tech_xmat: matrix of inputs for technology
# tech_ymat: matrix of outputs for technology
# eval_xmat: matrix of inputs for evaluating DMUs
# eval_ymat: matrix of outputs for evaluating DMUs
# Gx: direction of inputs
# Gy: direction of outputs
# convexity: convexity in the technology
# returns: returns to scale

ddf <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    Gx,
    Gy,
    convexity,
    returns
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to assess
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and number of outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  for (d in 1:eval_dmu) {

    # objective function
    objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
    objVal[1] <- 1

    # structure for lpSolve
    lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      add.constraint(lps, xt = c(Gx[, xi], tech_xmat[, xi]), "<=",  rhs = eval_xmat[d, xi])
    }

    # outputs
    for (yi in 1:nY) {
      add.constraint(lps, xt = c(- Gy[, yi], tech_ymat[, yi]), ">=", rhs =  eval_ymat[d, yi])
    }

    # technology
    if (returns == "variable") {

      if (convexity) {

        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)

      } else {

        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))

      }
    }

    set.bounds(lps, lower = c(- Inf, rep(0, tech_dmu)))

    solve(lps)
    scores[d, ] <- get.objective(lps)

  }

  return(scores)

}

# ================================= #
# Camano Dyson Luenberger Indicator #
# ================================= #

# tech_xmat_R: matrix of inputs for reference set to set technology.
# tech_ymat_R: matrix of outputs for reference set to set technology.
# eval_xmat_A_t: matrix of inputs for group A in "t" to evaluate.
# eval_ymat_A_t: matrix of outputs for group A in "t" to evaluate.
# eval_xmat_B_t: matrix of inputs for group B in "t" to evaluate.
# eval_ymat_B_t: matrix of outputs for group B in "t" to evaluate.
# Gx: direction of inputs.
# Gy: direction of outputs.
# convexity: convexity in the technology.
# returns: returns to scale.

# Formula (15)

cdli <- function (
    tech_xmat_R,
    tech_ymat_R,
    eval_xmat_A_t,
    eval_ymat_A_t,
    eval_xmat_B_t,
    eval_ymat_B_t,
    Gx,
    Gy,
    convexity,
    returns
  ) {

  # DDF: assess group B in "t" with reference set
  ddf_R_Bt <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_B_t,
    eval_ymat = eval_ymat_B_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in "t" with reference set
  ddf_R_At <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_A_t,
    eval_ymat = eval_ymat_A_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # ==
  # Camano Dyson Luenberger Indicator
  # ==

  CDLI <- mean(ddf_R_Bt) - mean(ddf_R_At)

  return (CDLI)

}

# ============== #
# Efficiency Gap #
# ============== #

# tech_xmat_A_t1: matrix of inputs for group A in t1 to set technology.
# tech_ymat_A_t1: matrix of outputs for group A in t1 to set technology.
# tech_xmat_B_t1: matrix of inputs for group B in t1 to set technology.
# tech_ymat_B_t1: matrix of outputs for group B in t1 to set technology.
# tech_xmat_A_t2: matrix of inputs for group A in t2 to set technology.
# tech_ymat_A_t2: matrix of outputs for group A in t2 to set technology.
# tech_xmat_B_t2: matrix of inputs for group B in t2 to set technology.
# tech_ymat_B_t2: matrix of outputs for group B in t2 to set technology.
# eval_xmat_A_t1: matrix of inputs for group A in t1 to evaluate.
# eval_ymat_A_t1: matrix of outputs for group A in t1 to evaluate.
# eval_xmat_B_t1: matrix of inputs for group B in t1 to evaluate.
# eval_ymat_B_t1: matrix of outputs for group B in t1 to evaluate.
# eval_xmat_A_t2: matrix of inputs for group A in t2 to evaluate.
# eval_ymat_A_t2: matrix of outputs for group A in t2 to evaluate.
# eval_xmat_B_t2: matrix of inputs for group B in t2 to evaluate.
# eval_ymat_B_t2: matrix of outputs for group B in t2 to evaluate.
# Gx: direction of inputs.
# Gy: direction of outputs.
# convexity: convexity in the technology.
# returns: returns to scale.

# Formula (19)

egap_ppli <- function (
    tech_xmat_A_t1,
    tech_ymat_A_t1,
    tech_xmat_B_t1,
    tech_ymat_B_t1,
    tech_xmat_A_t2,
    tech_ymat_A_t2,
    tech_xmat_B_t2,
    tech_ymat_B_t2,
    eval_xmat_A_t1,
    eval_ymat_A_t1,
    eval_xmat_B_t1,
    eval_ymat_B_t1,
    eval_xmat_A_t2,
    eval_ymat_A_t2,
    eval_xmat_B_t2,
    eval_ymat_B_t2,
    Gx,
    Gy,
    convexity,
    returns
    ) {

  # DDF: assess group B in t1 with technology B in t1
  ddf_Bt1_Bt1 <- ddf (
    tech_xmat = tech_xmat_B_t1,
    tech_ymat = tech_ymat_B_t1,
    eval_xmat = eval_xmat_B_t1,
    eval_ymat = eval_ymat_B_t1,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group B in t2 with technology B in t2
  ddf_Bt2_Bt2 <- ddf (
    tech_xmat = tech_xmat_B_t2,
    tech_ymat = tech_ymat_B_t2,
    eval_xmat = eval_xmat_B_t2,
    eval_ymat = eval_ymat_B_t2,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in t1 with technology A in t1
  ddf_At1_At1 <- ddf (
    tech_xmat = tech_xmat_A_t1,
    tech_ymat = tech_ymat_A_t1,
    eval_xmat = eval_xmat_A_t1,
    eval_ymat = eval_ymat_A_t1,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in t2 with technology A in t2
  ddf_At2_At2 <- ddf (
    tech_xmat = tech_xmat_A_t2,
    tech_ymat = tech_ymat_A_t2,
    eval_xmat = eval_xmat_A_t2,
    eval_ymat = eval_ymat_A_t2,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # ==
  # Efficiency Gap
  # ==

  EGAP <- mean(ddf_Bt1_Bt1) - mean(ddf_Bt2_Bt2) - (mean(ddf_At1_At1) - mean(ddf_At2_At2))

  return(EGAP)

}

# tech_xmat_A_t: matrix of inputs for group A in "t" to set technology.
# tech_ymat_A_t: matrix of outputs for group A in "t" to set technology.
# tech_xmat_B_t: matrix of inputs for group B in "t" to set technology.
# tech_ymat_B_t: matrix of outputs for group B in "t" to set technology.
# eval_xmat_A_t: matrix of inputs for group A in "t" to evaluate.
# eval_ymat_A_t: matrix of outputs for group A in "t" to evaluate.
# eval_xmat_B_t: matrix of inputs for group B in "t" to evaluate.
# eval_ymat_B_t: matrix of outputs for group B in "t" to evaluate.
# Gx: direction of inputs.
# Gy: direction of outputs.
# convexity: convexity in the technology.
# returns: returns to scale.

# Formula (16): EG

egap_cdli <- function (
    tech_xmat_A_t,
    tech_ymat_A_t,
    tech_xmat_B_t,
    tech_ymat_B_t,
    eval_xmat_A_t,
    eval_ymat_A_t,
    eval_xmat_B_t,
    eval_ymat_B_t,
    Gx,
    Gy,
    convexity,
    returns
    ) {

  # DDF: assess group B in "t" with technology B in "t"
  ddf_Bt_Bt <- ddf (
    tech_xmat = tech_xmat_B_t,
    tech_ymat = tech_ymat_B_t,
    eval_xmat = eval_xmat_B_t,
    eval_ymat = eval_ymat_B_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in "t" with technology A in "t"
  ddf_At_At <- ddf (
    tech_xmat = tech_xmat_A_t,
    tech_ymat = tech_ymat_A_t,
    eval_xmat = eval_xmat_A_t,
    eval_ymat = eval_ymat_A_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # ==
  # Efficiency Gap
  # ==

  EG <- mean(ddf_Bt_Bt) - mean(ddf_At_At)

  return (
    list (
      "EG" = EG,
      "ddf_A" = ddf_At_At,
      "ddf_B" = ddf_Bt_Bt
      )
  )
}

# ============== #
# Technology Gap #
# ============== #

# tech_xmat_R: matrix of inputs for reference set to set technology.
# tech_ymat_R: matrix of outputs for reference set to set technology.
# tech_xmat_A_t1: matrix of inputs for group A in t1 to set technology.
# tech_ymat_A_t1: matrix of outputs for group A in t1 to set technology.
# tech_xmat_B_t1: matrix of inputs for group B in t1 to set technology.
# tech_ymat_B_t1: matrix of outputs for group B in t1 to set technology.
# tech_xmat_A_t2: matrix of inputs for group A in t2 to set technology.
# tech_ymat_A_t2: matrix of outputs for group A in t2 to set technology.
# tech_xmat_B_t2: matrix of inputs for group B in t2 to set technology.
# tech_ymat_B_t2: matrix of outputs for group B in t2 to set technology.
# eval_xmat_A_t1: matrix of inputs for group A in t1 to evaluate.
# eval_ymat_A_t1: matrix of outputs for group A in t1 to evaluate.
# eval_xmat_B_t1: matrix of inputs for group B in t1 to evaluate.
# eval_ymat_B_t1: matrix of outputs for group B in t1 to evaluate.
# eval_xmat_A_t2: matrix of inputs for group A in t2 to evaluate.
# eval_ymat_A_t2: matrix of outputs for group A in t2 to evaluate.
# eval_xmat_B_t2: matrix of inputs for group B in t2 to evaluate.
# eval_ymat_B_t2: matrix of outputs for group B in t2 to evaluate.
# Gx: direction of inputs.
# Gy: direction of outputs.
# convexity: convexity in the technology.
# returns: returns to scale.

# Formula (20)

tgap_ppli <- function (
    tech_xmat_R,
    tech_ymat_R,
    tech_xmat_A_t1,
    tech_ymat_A_t1,
    tech_xmat_B_t1,
    tech_ymat_B_t1,
    tech_xmat_A_t2,
    tech_ymat_A_t2,
    tech_xmat_B_t2,
    tech_ymat_B_t2,
    eval_xmat_A_t1,
    eval_ymat_A_t1,
    eval_xmat_B_t1,
    eval_ymat_B_t1,
    eval_xmat_A_t2,
    eval_ymat_A_t2,
    eval_xmat_B_t2,
    eval_ymat_B_t2,
    Gx,
    Gy,
    convexity,
    returns
    ) {

  # DDF: assess group B in t1 with reference set
  ddf_R_Bt1 <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_B_t1,
    eval_ymat = eval_ymat_B_t1,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group B in t1 with technology B in t1
  ddf_Bt1_Bt1 <- ddf (
    tech_xmat = tech_xmat_B_t1,
    tech_ymat = tech_ymat_B_t1,
    eval_xmat = eval_xmat_B_t1,
    eval_ymat = eval_ymat_B_t1,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group B in t2 with reference set
  ddf_R_Bt2 <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_B_t2,
    eval_ymat = eval_ymat_B_t2,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group B in t2 with technology B in t2
  ddf_Bt2_Bt2 <- ddf (
    tech_xmat = tech_xmat_B_t2,
    tech_ymat = tech_ymat_B_t2,
    eval_xmat = eval_xmat_B_t2,
    eval_ymat = eval_ymat_B_t2,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in t1 with reference set
  ddf_R_At1 <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_A_t1,
    eval_ymat = eval_ymat_A_t1,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in t1 with technology A in t1
  ddf_At1_At1 <- ddf (
    tech_xmat = tech_xmat_A_t1,
    tech_ymat = tech_ymat_A_t1,
    eval_xmat = eval_xmat_A_t1,
    eval_ymat = eval_ymat_A_t1,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in t2 with reference set
  ddf_R_At2 <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_A_t2,
    eval_ymat = eval_ymat_A_t2,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in t2 with technology A in t2
  ddf_At2_At2 <- ddf (
    tech_xmat = tech_xmat_A_t2,
    tech_ymat = tech_ymat_A_t2,
    eval_xmat = eval_xmat_A_t2,
    eval_ymat = eval_ymat_A_t2,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # ==
  # Technology Gap
  # ==

  tgap <- (mean(ddf_R_Bt1 - ddf_Bt1_Bt1) + mean(ddf_Bt2_Bt2 - ddf_R_Bt2)) -
    (mean(ddf_R_At1 - ddf_At1_At1) + mean(ddf_At2_At2 - ddf_R_At2))

  return(tgap)

}

# tech_xmat_R: matrix of inputs for reference set to set technology.
# tech_ymat_R: matrix of outputs for reference set to set technology.
# tech_xmat_A_t: matrix of inputs for group A in year "t" to set technology
# tech_ymat_A_t: matrix of outputs for group A in year "t" to set technology
# tech_xmat_B_t: matrix of inputs for group B in year "t" to set technology
# tech_ymat_B_t: matrix of outputs for group B in year "t" to set technology
# eval_xmat_A_t: matrix of inputs for group A in year "t" to evaluate.
# eval_ymat_A_t: matrix of outputs for group A in year "t" to evaluate.
# eval_xmat_B_t: matrix of inputs for group B in year "t" to evaluate.
# eval_ymat_B_t: matrix of outputs for group B in year "t" to evaluate.
# Gx: direction of inputs.
# Gy: direction of outputs.
# convexity: convexity in the technology.
# returns: returns to scale.

# Formula (16: TG)

tgap_cdli <- function (
    tech_xmat_R,
    tech_ymat_R,
    tech_xmat_A_t,
    tech_ymat_A_t,
    tech_xmat_B_t,
    tech_ymat_B_t,
    eval_xmat_A_t,
    eval_ymat_A_t,
    eval_xmat_B_t,
    eval_ymat_B_t,
    Gx,
    Gy,
    convexity,
    returns
    ) {

  # DDF: assess group B in "t" with reference set
  ddf_R_Bt <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_B_t,
    eval_ymat = eval_ymat_B_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group B in "t" with technology B in "t"
  ddf_Bt_Bt <- ddf (
    tech_xmat = tech_xmat_B_t,
    tech_ymat = tech_ymat_B_t,
    eval_xmat = eval_xmat_B_t,
    eval_ymat = eval_ymat_B_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in "t" with reference set
  ddf_R_At <- ddf (
    tech_xmat = tech_xmat_R,
    tech_ymat = tech_ymat_R,
    eval_xmat = eval_xmat_A_t,
    eval_ymat = eval_ymat_A_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # DDF: assess group A in "t" with technology A in "t"
  ddf_At_At <- ddf (
    tech_xmat = tech_xmat_A_t,
    tech_ymat = tech_ymat_A_t,
    eval_xmat = eval_xmat_A_t,
    eval_ymat = eval_ymat_A_t,
    Gx = matrix(Gx, nrow = 1),
    Gy = matrix(Gy, nrow = 1),
    convexity = convexity,
    returns = returns
  )

  # ==
  # Technology Gap
  # ==

  tg <- mean(ddf_R_Bt) - mean(ddf_Bt_Bt) + mean(ddf_At_At) - mean(ddf_R_At)

  return (
    list (
      "tg" = tg,
      "DDF_A_R" = ddf_R_At,
      "DDF_A_A" = ddf_At_At,
      "DDF_B_R" = ddf_R_Bt,
      "DDF_B_B" = ddf_Bt_Bt
    )
  )
}