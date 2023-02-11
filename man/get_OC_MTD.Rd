\name{get_OC_MTD}
\alias{get_OC_MTD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate operating characteristics for finding the maximum tolerated dose (MTD)
}
\description{
Obtain the operating characteristics of the semi-mechanistic dose-finding design in oncology using pharmacokinetic/pharmacodynamic modeling }
\usage{
get_OC_MTD(ToxicityWindow = 28, TotaldoseMat = c(7,15,30,60,120), AdmTime = c(0, 7, 14, 21),
           Toxtarget = 0.3, CohortSize1 = 3, CohortSize2 = 3, MaxiSampleSize = 30, N1 = 30,
           espsilon1 = 0, espsilon2 = 0, psi1 = 0.1, psi2 = 0, start = 1000, thin = 1, update = 10000,
           SliceSampleSize = 2, mcmcsave = FALSE, mTPI = FALSE, IIV = TRUE, seed = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ToxicityWindow}{
time to observe whether patient has DLT or not
  }
  \item{TotaldoseMat}{
dosage at each level
  }
  \item{AdmTime}{
administered time for each dose level
  }
  \item{Toxtarget}{
desired toxicity level
  }
  \item{CohortSize1}{
cohortsize for phase 1 study
  }
  \item{CohortSize2}{
cohortsize for phase 2 study
  }
  \item{MaxiSampleSize}{
total sample size
  }
  \item{N1}{
sample size in phase 1 study
  }
  \item{espsilon1}{
a small value which makes (Toxtarget - espsilon1, Toxtarget + espsilon2) as an acceptable toxicity probability
  }
  \item{espsilon2}{
a small value which makes (Toxtarget - espsilon1, Toxtarget + espsilon2) as an acceptable toxicity probability
  }
  \item{psi1}{
early stopping criterion. If SafeProb[1]<=psi1 means no safe dose level
  }
  \item{psi2}{
admissable set criterion. If ifelse(SafeProb<psi2,Inf,ToxEst)
  }
  \item{start}{
burn in number for the posterior sampling
  }
  \item{thin}{
thinning number for the postrior sampling
  }
  \item{update}{
number of posterior sampling
  }
  \item{SliceSampleSize}{
number of slices in slice sampling
  }
  \item{mcmcsave}{
whether to save mcmc result
  }
  \item{mTPI}{
whether to use mTPI method
  }
  \item{IIV}{
whether the population variables come from a prior distribution
  }
  \item{seed}{
seed for rerun the code
  }


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Su, X., Li, Y., Mueller, P., Hsu, C. W., Pan, H. & Do, K. A. A semi-mechanistic dose-finding design in oncology using pharmacokinetic/pharmacodynamic modeling. Pharmaceutical Statistics. 2022 Jun. doi.org/10.1002/pst.2249
}
\author{
Chia-Wei Hsu, Xiao Su, Haitao Pan, Yisheng Li
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
get_OC_MTD(ToxicityWindow = 28, TotaldoseMat = c(7,15,30,60,120), AdmTime = c(0, 7, 14, 21),
           Toxtarget = 0.3, CohortSize1 = 3, CohortSize2 = 3, MaxiSampleSize = 30, N1 = 30,
           espsilon1 = 0, espsilon2 = 0, psi1 = 0.1, psi2 = 0, start = 1000, thin = 1, update = 10000,
           SliceSampleSize = 2, mcmcsave = FALSE, mTPI = FALSE, IIV = TRUE, seed = 6)
}

