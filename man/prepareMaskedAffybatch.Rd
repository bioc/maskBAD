% --- Source file: man/prepare.CDF.fromTable.Rd ---
\name{prepareMaskedAffybatch}
\alias{prepareMaskedAffybatch}
\title{ Creating a new CDF }
\description{ Create a new affyBatch, withprobes and probesetsdefined by
  mask.
}
\usage{
prepareMaskedAffybatch(affy,cdfTablePath,exmask="none",cdfName="new_cdf",exclude=NA,cutoff=0.2)
}
\arguments{
  \item{affy}{ An object of class \code{AffyBatch}.}
  \item{cdfTablePath}{ Location of the probe information table. This is
    a plain text file with probes to build new cdf. It should contain 3
    or 5 columns. Column 1: Probeset ID. Column 2: probe
    x-coordinate. Column 3: probes y-coordinate. Optional column 4:
    Mismatch probe x-coordinate. Optional column 5: Mismatch probe y
    coordinate. }
  \item{exmask}{ Data frame with probe information, for example
    first element of the output of function \code{mask}. Should contain:
    column 1: probe x-coordinate, column 2:probe y coordinate, column 3
    :probeset, column 4: quality score: values to based filtering on,
    probes with values smaller than cutoff are discarded. }
  \item{cdfName}{ Name for the new CDF. }
  \item{cutoff}{ With \code{mask.object}, defines the minimum quality
    score necessary for a probe to qualify to the new cdf. }
  \item{exclude}{ Default 'NA'. If exclude set to a number>0, probesets
    with less than 'exclude' probes remaining after masking are excluded
    from the new affyBatch object. }
}

\details{ The function \code{ prepareMaskedAffybatch } creates a new
  affyBatch including only the probes remaining after masking. Set of
  probes might be defined by a txt file, with \code{cdfTablePath}
  argument, or by a data frame \code{mask.object} and \code{cutoff} the
  probes have to exceed to be used in the new cdf.
}

\value{
  \item{ newAffyBatch }{ A list with an affyBatch object and an
    environment for the new CDF identifier. }
}
\references{ Dannemann et al, {The effects of probe binding
    affinity differences on gene expression measurements and how to deal
    with them}. Bioinformatics 2009 } %% to do
\author{ Michael Lachmann, Mehmet Somel, Michael Dannemann, Anna Lorenc }
\seealso{ \code{\link{mask}, \link{overlapExprExtMasks}, \link{plotProbe}}}
\examples{
## prepare new affy batch after masking
## using the expression mask object from the example of the mask function
data(AffyBatch)
data(exmask)
## AffyBatch object before masking
newAffyBatch
affyBatchAfterMasking <-
   prepareMaskedAffybatch(affy=newAffyBatch,exmask=exmask$probes)
## AffyBatch object after masking
affyBatchAfterMasking
}
\keyword{internal}
