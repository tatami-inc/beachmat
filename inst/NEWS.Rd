\name{beachmatnews}
\title{beachmat News}
\encoding{UTF-8}

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

