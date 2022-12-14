% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterAnnotation.r
\name{clusterAnnotation}
\alias{clusterAnnotation}
\title{Cluster annotation}
\usage{
clusterAnnotation(file, attribute, clus, outputDir)
}
\arguments{
\item{file}{meta data file}

\item{attribute}{Field one from meta data}

\item{clus}{Field two from meta data}

\item{outputDir}{output directory for cluster annotation file}
}
\value{
PDF file with overlap between two annotation fields
}
\description{
This function allows the users to compare different annotation classes to see the overlap between two meta fields
}
