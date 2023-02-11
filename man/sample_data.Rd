\name{sample_data}
\alias{sample_data}
\docType{data}
\title{
Sample dataset
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(sample_data)}
\format{
  The sample_data list contains two elements. The first element is the toxicity data frame (30 rows and 9 columns) by each subject at the end of follow-up time (28 days)
  The second element is the PK data frame (180 rows and 10 columns) of 6 measurements of drug concentration in the human body wihtin a day at different time points for each patient after the first dosage.
  The variables are as follows:
  \describe{
   \item{subject}{The number identifying the subject on whom the observation is made.}
  	\item{visit}{The number indicating the order of the sample for each patient.}
  	\item{meanC}{The mean concentration of the dosage in the human body}
  	\item{logmeanC}{The logarithm of the meanC}
  	\item{epsilonC}{Measurement error of the logarithm of the meanC}
  	\item{Cobs}{Observed concentration of the dosage in the human body}
  	\item{logCobs}{The logarithm of the observed concentration of the dosage in the human body}
  	\item{t}{Measurement time. The unit is day}
  	\item{Cnobs}{Number of observation for the subject}
  	\item{regime}{Dose level}
  	\item{delta}{Whether you have DLT or not}
  	\item{ToxRate}{The probability to get toxicity}
  	\item{ToxicityWindow}{The day to obsserve whether you have DLT or not}
  	\item{k}{The elimination rate of the person}
  	\item{v}{The distribution volume of the person}
  	\item{logk}{The logarithm of the k}
  	\item{logv}{The logarithm of the v}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Su, X., Li, Y., Mueller, P., Hsu, C. W., Pan, H. & Do, K. A. A semi-mechanistic dose-finding design in oncology using pharmacokinetic/pharmacodynamic modeling. Pharmaceutical Statistics. 2022 Jun. doi.org/10.1002/pst.2249
}

\seealso{
  \code{\link{next_MTD}},
  \code{\link{select_MTD}},
  \code{\link{get_OC_MTD}},
}



\examples{
data(sample_data)
head(sample_data)
}
\keyword{datasets}
