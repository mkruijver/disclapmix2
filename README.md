# disclapmix2

An experimental R-package for estimation of forensic Y-STR haplotype frequencies that generalises the model implemented in the [disclapmix](https://github.com/mikldk/disclapmix) package. This extension allows for multi-copy loci, partial repeats and null alleles. For integer valued data, the results should be concordant with the original `disclapmix` implementation. Unlike the original implementation, this package uses an off-the-shelf solver for numerical likelihood optimisation.
