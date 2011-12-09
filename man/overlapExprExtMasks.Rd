% --- Source file: man/errors.Rd ---
\name{overlapExprExtMasks}
\alias{overlapExprExtMasks}
\title{ Error Analysis of Masking Results }
\description{ Expression mask results for a range of cutoff values are
  compared with an external  mask (for example a mask based on sequence
  data) and  type 1 and type 2 errors are estimated. }
\usage{
overlapExprExtMasks(probes,seqdata,cutoffs="none",wilcox.ks=FALSE,sample=10,plotCutoffs=TRUE,verbose=TRUE)
}
\arguments{
  \item{probes}{ A matrix with 3 columns. The first and second column
    represent the x and y coordinates on the Microarray. The third column contains
    a quality entry for each probe, e.g. the quality score obtained
    from \code{ mask } analysis.}
  \item{seqdata}{ A matrix with 3 columns containing x, y coordinates and
    0,1 entries in column 3, defining whether a probe
    has a sequence difference (0) or not. }
  \item{cutoffs}{ A vector including all cutoff values for the quality
    scores of an expression mask that should be used for the error
    analysis. If no cutoffs are given (default is "none") the cutoffs
    are the quantiles of the quality scores starting from 0 to 1 in
    steps of 0.01. }
  \item{wilcox.ks}{Logical, default=FALSE element determining whether the Kolmogorow-Smirnow
    Test and Wilcoxon Rank Test analysis should be performed (see reference below).}
  \item{sample}{ To compare the p value distribution with the
    Kolmogorow-Smirnow Test and Wilcoxon Rank Test for
    different cutoffs, the sampling option can be used to compute the quality
    score distribution for different cutoffs. This value indicates how often
    the sampling should be performed. }
  \item{plotCutoffs}{Logical, default=TRUE element determining whether the cutoffs
    should be drawn in the overlap plot.}
  \item{verbose}{ Logical. If 'TRUE', it writes out some messages
    indicating progress. If 'FALSE' nothing should be printed.}
}
\details{
  The function \code{ overlapExprExtMasks } compares expression mask
  results with an external (for example sequence-based) mask and might
  help to choose a quality score cutoff for masking probes.
}
\value{
A list of five objects will be returned.
  \item{ type1 }{ A vector of the type 1 error for each cutoff. }
  \item{ type2 }{ A vector of the type 2 error for each cutoff. }
  \item{ confT1 }{A matrix with the upper (column 1) and lower (column
    2) confidence intervals for the type 1 error.}
  \item{ confT2 }{A matrix with the upper (column 1) and lower (column
    2) confidence intervals for the type 2 error.}
  \item{ ksP }{  If wilcox.ks is 'TRUE', a vector of quality scores
    from a two sample Kolmogorov-Smirnov comparing distributions 
    of quality score for probes designated as BAD and not in external
    mask. } 
  \item{ wilcoxonP }{  If wilcox.ks is 'TRUE', a vector of quality
    scores from a two sample Wilcoxon rank test comparing distributions
    of quality score for probes designated as BAD and not in external
    mask. }
  \item{ ksBoot }{ For each cutoff sample(default=10) times cutoff
    values for the Kolmogorov-Smirnov test will be generated.}
  \item{ wilcoxBoot }{ For each cutoff sample(default=10) times cutoff
    values for the wilcoxon rank sum test will be generated.}
  \item{cutoffs}{List of cutoffs used for the error analysis}
  \item{testCutoffs}{ If wilcox.ks is 'TRUE', a list with cutoff
    information will be provided. The first list entry includes all
    cutoffs used in the two sample Kolmogorov-Smirnov test and the two
    sample wilcoxon rank sum test analysis will be produced. A cutoff
    can appear sample(default=10) times. In theory there should be
    sample times the number of cutoff values entries in this vector, but
    usually there are fewer entries, because for certain cutoff values,
    it is not possible to calculate the exact p value in one of the
    tests. The second list entry transforms the cutoffs in ranks and can
    be used for the plotting of the test results. }
}
\references{ Dannemann et al, {The effects of probe binding
    affinity differences on gene expression measurements and how to deal
    with them}. Bioinformatics 2009 } %% to do
\author{ Michael Dannemann }
\seealso{ \code{\link{mask}, \link{prepareMaskedAffybatch},
    \link{plotProbe} } }
\examples{
## loading mask on all genes (exmask1) of the same dataset
data(exmask)
overlapExSeq <- overlapExprExtMasks(exmask$probes[,1:3],sequenceMask[,c(1,2,4)])

## plot results
plot(overlapExSeq$type1,overlapExSeq$type2,type="l",col="red",
     main="Overlap expression based mask - sequence based mask",xlab="Type 1",ylab="Type 2")
abline(1,-1,col="gray")

## performing wilcoxon rank sum test and Kolmogorov-Smirnov test on
## expression mask with all genes (exmask)
overlapTests <-
  overlapExprExtMasks(exmask$probes[,1:3],sequenceMask[,c(1,2,4)],wilcox.ks=TRUE)
layout(matrix(1:2,ncol=1))
plot(overlapTests$testCutoff[[1]],overlapTests$ksBoot,col="red",main="Kolmogorov-Smirnov Test",xlab="Quality score cutoff",
     ylab="p value (Kolmogorov-Smirnov Test)",ylim=c(0,1),pch=16,xaxt="n")
axis(1,at=1:length(unique(overlapTests$testCutoff[[2]])),labels=signif(unique(overlapTests$testCutoff[[2]]),2),las=3)
lines(which(unique(overlapTests$testCutoff[[2]]) \%in\% overlapTests$testCutoff[[2]]),overlapTests$ksP[!is.na(overlapTests$ksP)],type="p",pch=16,cex=0.8)
plot(overlapTests$testCutoff[[1]],overlapTests$wilcoxonBoot,col="green",main="Wilcoxon Rank Sum Test",xlab="Quality score cutoff",
     ylab="p value (Wilcoxon Rank Sum Test)",ylim=c(0,1),pch=16,xaxt="n")
axis(1,at=1:length(unique(overlapTests$testCutoff[[2]])),labels=signif(unique(overlapTests$testCutoff[[2]]),2),las=3)
lines(which(unique(overlapTests$testCutoff[[2]]) \%in\% overlapTests$testCutoff[[2]]),overlapTests$wilcoxonP[!is.na(overlapTests$wilcoxonP)],type="p",pch=16,cex=0.8)

}
\keyword{internal}
