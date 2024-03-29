% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_haplotypes.R
\name{haplotype_counts}
\alias{haplotype_counts}
\title{Count the number of times each haplotype occurs}
\usage{
haplotype_counts(x)
}
\arguments{
\item{x}{DataFrame (by locus) of character vectors containing haplotypes (rows) where alleles are separated by comma's, e.g. "13,14.2" is a haplotype}
}
\value{
Integer vector with count for each row in DataFrame
}
\description{
Count the number of times each haplotype occurs
}
\examples{
# read haplotypes
h <- readxl::read_excel(system.file("extdata","South_Australia.xlsx",
package = "disclapmix2"), 
col_types = "text")[-c(1,2)]

# obtain counts
counts <- disclapmix2::haplotype_counts(h)

# all haplotypes in the dataset are unique
stopifnot(all(counts == 1))
}
