\name{symbolize}
\alias{symbolize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert original data to symbols}
\description{An auxiliary function of NTest that converts original data into symbolic data.}
\usage{symbz <- symbolize(x,y,z,m)}
\arguments{
  \item{x}{data, dependent at t}
  \item{y}{data, dependent at t-1}
  \item{z}{data, independent}
  \item{m}{lag scope}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{symbz = a structure with sorted indices of mhx, mhy, mhz, mhxy, mhxz and mhxyz}
\references{}
\author{
Andres Romeu
    Fundamentos del Analisis Economico
    Universidad de Murcia 30100 SPAIN}
\note{Translated from original Matlab code by my colleague Ruiz-Marin}

\seealso{TEntropy}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
