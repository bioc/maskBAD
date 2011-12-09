% --- Source file: man/mask.Rd ---
\name{mask}
\alias{mask}
\title{ Filtering/Masking expression data }
\description{ Identyfying probes with binding affinity difference (BAD probes) between two groups of samples on the basis of expression data. }
\usage{
mask(affy,exprlist=NULL,useExpr=TRUE,ind,PM=FALSE,verbose=TRUE)
}
\arguments{
  \item{affy}{An object of class \code{AffyBatch}.}
  \item{exprlist}{ A vector with  probesetnames
    to be used. If NULL, all probesets are analyzed. }
  \item{useExpr}{ Logical. If  'TRUE', only expressed genes (see
    Details) are used. If  'FALSE', all probes are analyzed. }
  \item{ind}{ Numeric vector, with values 1 and 2, definining group
    assignement for samples in \code{affy}. }
  \item{PM}{ Logical. If 'TRUE', only probes with a mean pm
    value greater than the mean mm value are used.}
  \item{verbose}{ Logical. If 'TRUE', it writes out some messages
    indicating progress. If 'FALSE' nothing should be printed.}
}
\details{
  The function \code{mask} identifies in expression data probes which
  binding affinity (BAD  probes) differs between two groups of samples,
  e.g two species. The basic input data is \code{AffyBatch
    object}(expresssion data prepared using the function \code{ ReadAffy
  } from the library \code{ Affy}) and a vector defining group
  assignement of samples. As masking is based on expression values, only
  expressed probes should be used. As a default they are defined by the
  \code{ affy } function \code{ mas5calls} and condition of being
  expressed (having "P" value) in at least 90\% of samples from each
  group, but any set of probesets might be submitted with \code{
    exprlist} argument.  Probes are analyzed for difference in binding
  affinity between groups. Each probe is assigned a quality score, based
  on  all pairwise probes' correlations within probesets (for details
  see vignette or paper).  Probes' quality scores, their x and y
  coordinates on the microarray and the probeset  names are stored in a
  matrix. }
\value{
A list of two objects will be returned.
  \item{ probes }{A data frame with x,y coordinates, quality score and probeset for each analyzed probe.}
  \item{ notUsed }{If PM=TRUE: A vector with unused probes having a lower pm
    mean value than mm mean value.}
}
\references{ Dannemann et al, {The effects of probe binding
    affinity differences on gene expression measurements and how to deal
    with them}. Bioinformatics 2009 \\
  Khaitovich et al,Parallel Patterns of Evolution in the
  Genomes and Transcriptomes of Humans and Chimpanzees, Science 2005} %% to do
\author{ Michael Dannemann, Michael Lachmann }
\seealso{ \code{\link{overlapExprExtMasks}, \link{prepareMaskedAffybatch},
    \link{mas5calls}, \link{plotProbe} } }
\examples{
data(AffyBatch)
## we provide 20 samples (10 for both human and chimpanzee)
## the first 10 entries are chimpanzee samples the last 10 from human
ind.vec=rep(1:2,each=10)
## mask on AffyBatch with all genes 
exmask <-
mask(newAffyBatch,ind=ind.vec,PM=TRUE,useExpr=FALSE)
}

\keyword{internal}
