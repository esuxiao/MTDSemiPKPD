\name{next_MTD}
\alias{next_MTD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a maximum tolerated dose (MTD)
}
\description{
Determine the dose for the next cohort of new patients for single-agent trials that aim to find a MTD
}
\usage{
next_MTD(data = data, ToxicityWindow = 28, AdmTime = c(0, 7, 14, 21), TotaldoseMat = c(7, 15, 30, 60, 120),
         Toxtarget = 0.3, espsilon1 = 0, espsilon2 = 0, psi1 = 0.5, psi2 = 0.6, mTPI = FALSE, IIV = TRUE,
         update = 10000, start = 1000, thin = 1, mcmcsave = FALSE, seed = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
data to be used for the analysis
}
 \item{ToxicityWindow}{
time to observe whether patient has DLT or not
}

 \item{AdmTime}{
administered time for each dose level
}
 \item{TotaldoseMat}{
dosage at each level
}
 \item{Toxtarget}{
desired toxicity level
}
 \item{espsilon1}{
a small value which makes (Toxtarget - espsilon1, Toxtarget + espsilon2) as an acceptable toxicity probability
}
 \item{espsilon2}{
a small value which makes (Toxtarget - espsilon1, Toxtarget + espsilon2) as an acceptable toxicity probability
}
 \item{psi1}{
early stopping criterion
}
 \item{psi2}{
admissable set criterion
}

 \item{mTPI}{
whether to use mTPI method
}

 \item{IIV}{
whether the population variables come from a prior distribution
}

 \item{update}{
number of posterior sampling
}

\item{start}{
burn in number for the posterior sampling
}

\item{thin}{
thinning number for the postrior sampling
}

\item{mcmcsave}{
Whether save your MCMC result
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
data(sample_data)
next_MTD(data = sample_data, ToxicityWindow = 28, AdmTime = c(0, 7, 14, 21), TotaldoseMat = c(7, 15, 30, 60, 120),
         Toxtarget = 0.3, espsilon1 = 0, espsilon2 = 0, psi1 = 0.5, psi2 = 0.6, mTPI = FALSE, IIV = TRUE,
         update = 10000, start = 1000, thin = 1, mcmcsave = FALSE, seed = 10)
}

