% --- Source file: man/seq.mask.Rd ---
\name{plotProbe}
\alias{plotProbe}
\title{ Plot probes }
\description{ Pairwise plot probes of a probeset.
}
\usage{
plotProbe(affy,probeset,probe=NA,probeXY=NA,scan=TRUE,ind,exmask="none",seqmask="none",names=FALSE)
}
\arguments{
  \item{affy}{ An object of class \code{AffyBatch}. }
  \item{probeset}{Probe set name (Affymetrix ID).}
  \item{probe}{ Number of the main probe. }
  \item{probeXY}{ If probe is NA the x and y coordinates of the main
    probe can be given in the format 'x.y'. }
  \item{scan}{ If scan is 'TRUE', each probewise comparison of the probe
    against all other probes in this probeset will be performed
    seperately. If scan is 'FALSE', all plots will be plotted in one
    layout. The layout has 3 columns. If the number of remaining probes
    that the probe should be compared with is not a multiple of 3, the
    number of probes will be reduced to the next lower multiple of 3. }
  \item{ind}{ Numeric vector, with values 1 and 2, definining group
    assignement for samples in \code{affy}.  }
  \item{exmask}{ Optional: an expression mask object for this affy
    batch. Data frame with probe information, for example first element
    of the output of function \code{mask}. Should contain: column 1:
    probe x-coordinate, column 2:probe y coordinate, column 3: probeset,
    column 4: quality score: values to based filtering on, probes with
    values smaller than cutoff are discarded.
 }
  \item{seqmask}{ Optional: a sequence mask object for this
    mask. }
  \item{names}{ If 'TRUE' , the sample names
    are plotted to identify each individual. }
}

\details{ The function \code{plotProbe} plots single probe against all
  other probes of its probe set. The information from the expression
  based mask, the sequence based mask and the test for the two plotted
  probes is shown.
}

\references{ Dannemann et al, {The effects of probe binding
    affinity differences on gene expression measurements and how to deal
    with them}. Bioinformatics 2009 } %% to do
\author{ Michael Dannemann }
\seealso{ \code{\link{mask}, \link{overlapExprExtMasks},
    \link{prepareMaskedAffybatch}}}
\examples{
data(exmask)
data(AffyBatch)
## plot for one probe comparisons with other probes of the probeset
## for a random probeset
availableProbesets <- as.character(unique(exmask$probes[,4]))
availableProbesets
## scan the plots
\dontrun{plotProbe(affy=newAffyBatch,probeset=availableProbesets[22],probe=5,scan=TRUE,ind=rep(1:2,each=10),exmask=exmask$probes,seqmask=sequenceMask)}
## scan with names=TRUE
\dontrun{plotProbe(affy=newAffyBatch,probeset=availableProbesets[22],probe=5,scan=TRUE,ind=rep(1:2,each=10),exmask=exmask$probes,seqmask=sequenceMask,names=TRUE)}
## plot with given x y information
\dontrun{plotProbe(affy=newAffyBatch,probeset=availableProbesets[22],probeXY="313.415",scan=TRUE,ind=rep(1:2,each=10),exmask=exmask$probes,seqmask=sequenceMask,names=TRUE)}
## all plots in one layout
plotProbe(affy=newAffyBatch,probeset=availableProbesets[22],probe=5,scan=FALSE,ind=rep(1:2,each=10),exmask=exmask$probes,seqmask=sequenceMask)}

\keyword{internal}
