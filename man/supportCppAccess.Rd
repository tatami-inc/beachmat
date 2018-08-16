\name{supportCppAccess}
\alias{supportCppAccess}
\alias{supportCppAccess,ANY-method}

\title{Support C++ access}
\description{Does the current matrix class support C++ access?}

\usage{
\S4method{supportCppAccess}{ANY,character}(x)
}

\arguments{
\item{x}{A matrix-like object.}
}

\value{
A logical scalar specifying whether C++ access is supported for \code{x}.
}

\details{
This function is called by the \pkg{beachmat} C++ API upon encountering an unknown matrix type.
If it returns \code{TRUE}, we assume that the package used to define the class of \code{x} also contains registered C++ functions to access rows or columns of \code{x}.
This allows \pkg{beachmat} to use those C++ functions to directly access values of \code{x}, rather than relying on block realization of an unknown matrix.

If you need to use this function, you almost certainly need to read the vignette on \dQuote{Extending beachmat}.
}

\author{
Aaron Lun
}

\examples{
library(Matrix)
Y <- Matrix(0, 10, 10)    
supportCppAccess(Y)

# Note that certain matrix types are always supported
# by beachmat but still return FALSE.
X <- matrix(0, 10, 10)
supportCppAccess(X)
}
