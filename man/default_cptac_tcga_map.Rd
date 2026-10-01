\name{default_cptac_tcga_map}
\alias{default_cptac_tcga_map}
\title{Default CPTAC to TCGA Cancer Mapping}
\usage{
default_cptac_tcga_map() 

}
\description{
Returns the manually reviewed mapping between the curated CPTAC cancer labels
used by the proteogenomic validation workflow and TCGA cancer abbreviations.
}
\details{
This mapping is used only for cancer-level and effect-size-level comparison.
It does not imply patient-level matching between TCGA and CPTAC.
}
\value{
A data.table with columns \code{cptac_cancer}, \code{tcga_cancer}, and
\code{note}.
}
\examples{
default_cptac_tcga_map()
}
