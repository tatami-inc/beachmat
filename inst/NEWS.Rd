\name{beachmatnews}
\title{beachmat News}
\encoding{UTF-8}

\section{Version 2.24.0}{\itemize{
\item Automatically cast integer/logical NAs to their double equivalents in the matrices wrapped by \code{initializeCpp()}.
This can be disabled by setting \code{.check.na=FALSE}.

\item Added a \code{initializeCpp()} method for the ConstantArraySeed class. 

\item Using a \code{TileDBArraySeed} in \code{initializeCpp()} now automatically dispatches to \pkg{beachmat.tiledb}.
}}

\section{Version 2.22.0}{\itemize{
\item Moved all C++ libraries to the \pkg{assorthead} package.
This uses the latest versions of the \pkg{tatami} framework.

\item Automatically attempt to use \pkg{beachmat.hdf5} (if installed) when encountering the various \pkg{HDF5Array} classes.

\item Minor optimization when dealing with some known no-op matrices, e.g., \code{WrapperArraySeed} objects from \pkg{alabaster.matrix}.

\item Map delayed type coercions from \code{type<-} to their corresponding native representations (or no-ops).

\item Bugfix for \code{initializeCpp()} to work properly with dense \pkg{Matrix} instances.

\item Soft-deprecated \code{whichNonZero()} in favor of \pkg{SparseArray}'s new \code{nzwhich()} and \code{nzvals()} functions.

\item Exported \code{tatami.*} utilities for manipulating already-initialized pointers from \code{initializeCpp()}.
This includes all of the delayed operation wrappers, some previously internal functions for extracting matrix data,
a new function for realizing the matrix into a standard R representation,
and a new function for matrix multiplication.
}}

\section{Version 2.14.0}{\itemize{
\item Added the \pkg{tatami} C++ library for LinkingTo from other packages.
}}

\section{Version 2.8.0}{\itemize{
\item Improve the efficiency of sparse row subsetting in non-DelayedArray \code{rowBlockApply()}.

\item Avoid overhead of \pkg{DelayedArray} block processing when DelayedArray is pristine and the type is supported.

\item Migrated \code{whichNonZero()} from \pkg{scuttle}.

\item Added \code{toCsparse()} to make it easier to convert SparseArraySeeds to CsparseMatrixes.

\item Added \code{realizeFileBackedMatrix()} to, well, realize a DelayedMatrix with file-backed components.
}}

\section{Version 2.6.0}{\itemize{
\item Added version 3 of the C++ API.
This relies on \pkg{DelayedArray} block processing to handle the various abstractions,
and provides convenient C++ abstractions for the different input types.
More methods for non-zero extraction are provided for sparse algorithms,
and the class hierarchy handles flexibility of the input types rather than requiring user-side templating.

\item Added the \code{colBlockApply()} and \code{rowBlockApply()} functions for more convenient and efficient column/row-block processing.

\item Added doxygen documentation for the C++ API.
}}

\section{Version 2.0.0}{\itemize{
\item Converted beachmat into a header-only library to simplify linking, avoid ABI incompatibilities.

\item Minor API change for external native support to avoid using C++ classes with C-style linkage.
}}

\section{Version 1.4.0}{\itemize{
\item Removed native support for RleMatrix and packed symmetric matrices.

\item Added multi-row/column getters.

\item Added mechanism for native support of arbitrary developer-defined matrices.

\item Switched to row/colGrid() for defining chunks in unsupported matrices.
}}

\section{Version 1.2.0}{\itemize{
\item Changed environment variable to BEACHMAT_RPATH for consistency with other packages.

\item Added native support for transposition and subsetting in DelayedMatrix objects.

\item Added support for chunk-by-chunk realization of otherwise unsupported matrices, including DelayedMatrix objects with other delayed operations.

\item Added the get_const_col_indexed() method for input matrices, especially fast for sparse representations. 

\item Added the set_col_indexed() and set_row_indexed() methods for output matrices.

\item Updated vignettes and expanded the test suite.
}}

\section{Version 1.0.0}{\itemize{
\item
New package beachmat, which provides a C++ API for handling various R matrix types.
}}

