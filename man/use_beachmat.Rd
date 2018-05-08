\name{use_beachmat}
\alias{use_beachmat}
\title{Use beachmat}
\usage{
use_beachmat()
}
\description{
Creates \code{src/}, adds needed packages to \code{DESCRIPTION}, and prompts user to create or edit files required for linking to \pkg{beachmat} with C++ code from another package.
}
\details{
For more details see the 'Linking' vignette, available by running \code{vignette("linking", package = "beachmat")}.
}
\value{
Called for its side effect of creating files to set up a package for use with \pkg{beachmat}.
Assumes that the current working directory is the root directory of the package under development.
}

\author{
Peter Hickey
}
\examples{
\dontrun{
beachmat:::use_beachmat()
}
}
