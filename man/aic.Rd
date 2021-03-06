% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{aic}
\alias{aic}
\title{Akaike Information Criterion}
\usage{
aic(table, dimension, clusters, scheme)
}
\arguments{
\item{table}{object of class "table".}

\item{dimension}{integer. Whether to use rows (1) or columns (2). Default is 1.}

\item{clusters}{list of integer vectors. Each vector should define a cluster
by specifing row or column indices of its memebrs. Clusters must not overalap.}

\item{scheme}{character. Name of the assumed table's generation scheme.
Can be either "multinomial" or "independent multinomial"}
}
\value{
double. Value of AIC statisctic.
}
\description{
See section 4.1 in thesis, for details
}
\examples{
clusters <- clusters2list(stats::cutree(hclust_table(ksarakil, 1), k = 3))
aic(ksarakil, 1, clusters, "independent multinomial")
}
