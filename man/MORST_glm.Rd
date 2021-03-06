% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MORST_glm.R
\name{MORST_glm}
\alias{MORST_glm}
\title{The score version of MORST for Generalized Linear Models.}
\usage{
MORST_glm(
  X,
  obj,
  alpha = 0.05,
  weights = NULL,
  tau.type = "approx",
  target_power = 0.5,
  n.points = 50
)
}
\arguments{
\item{X}{a numeric matrix or dgCMatrix of predictors.}

\item{obj}{an output from \code{Null_model_glm}.}

\item{alpha}{the alpha parameter in MORST. It is suggested to be the significance level.}

\item{weights}{a numeric vector of nonnegative weights. If \env{NULL}, the equal weight is used.}

\item{tau.type}{either "minimax" or "approx". See \code{\link{tau_c}} for details.}

\item{target_power}{a value that is used when \emph{tau.type == "approx"}. See \code{\link{tau_c}} for details.}

\item{n.points}{number of grid points used when \emph{tau.type == "minimax"}. See \code{\link{tau_c}} for details.}
}
\value{
The p-value of MORST.
}
\description{
Calculate the score version MORST p-value for GLM. For genetic association study, please use the function \code{\link{SetBasedTests}}.
}
\examples{
X<-matrix(rnorm(20000),ncol=20); Z=matrix(rnorm(nrow(X)*4),ncol=4)

### linear regression for continuous outcome
Y<-rnorm(nrow(X));obj<-Null_model_glm(Y,Z,family="gaussian")
MORST_glm(X,obj)

### Logistic regression for binary outcome
Y<-rbinom(nrow(X),1,0.4);obj<-Null_model_glm(Y,Z,family="binomial")
MORST_glm(X,obj,alpha = 1e-04)

### Binomial outcome
Y<-rbinom(nrow(X),5,0.4);Y<-cbind(Y,5-Y);obj<-Null_model_glm(Y,Z,family="binomial")
MORST_glm(X,obj,weights=runif(ncol(X)))

### Poisson outcome
Y<-rpois(nrow(X),5);obj<-Null_model_glm(Y,Z,family="poisson")
MORST_glm(X,obj,alpha = 1e-04,weights=runif(ncol(X)),tau.type = "minimax")
}
\references{
Liu, Y., Li, Z., and Lin, X. (2020+) A Minimax Optimal Ridge-Type Set Test
for Global Hypothesis with Applications in Whole Genome Sequencing Association Studies.
\emph{Journal of the American Statistical Association}. Accepted.
}
\author{
Yaowu Liu
}
