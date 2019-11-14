\name{NTest}
\alias{NTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Net Transfer Entropy Causality Test for Panel Data}
\description{Computes the Net Transfer Entropy symbolic test in Camacho, Romeu and Ruiz-Marin (2019)}
\usage{nte <- NTest(data,noms,m,r,B)}
\arguments{
  \item{data}{pdataframe, first column Y second column X}
  \item{noms}{Vector of names of variables to test}
  \item{m}{Length of window}
  \item{r}{Lags to include}
  \item{B}{Number of bootstraps}
}
\details{}
\value{nte = a structure with
  \item{testXcY}{test statistic of X causes Y}
  \item{testYcX}{test statistic of Y causes X}
  \item{NETransfer}{net transfer statistic. Positive: XcY > YcX. Negative: YcX > XcY}
  \item{BTS}{a data frame with boots of XcYtest,YcXtest and NTEtest}
}
\references{Camacho, Romeu and Ruiz-Marin (2019)}
\author{Andres Romeu
    Fundamentos del Analisis Economico
    Universidad de Murcia 30100 SPAIN
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
