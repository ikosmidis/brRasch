% -*- coding: utf-8 -*-
\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options: width=0.5in}}}
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\sspace}{\ifelse{latex}{\out{~}}{ }}
%% some pre-defined commands: \R, \code, \acronym, \url, \file, \pkg

\name{NEWS}
\title{News for package \CRANpkg{brRasch}}
\encoding{UTF-8}

\section{Changes in brRasch version 0.1.1.9000 (2015-12-28)}{
  \subsection{ENHANCEMENTS}{
    \itemize{
      \item Improvements in the plotting functions
      \item 1PL models can now be fit by setting dim = 0
      \item Improvements in the plotting functions (use of ggplot2
      throughout)
      \item Enhancements in the print.brRasch function
    }
  }
}


\section{Changes in brRasch version 0.1.0.9000 (2015-03-29)}{
  \subsection{ENHANCEMENTS}{
    \itemize{
      \item Added the option for reguralized estimation and implemented
      L2 penalty on the discrimination parameters.
    }
  }
}


\section{Changes in brRasch version 0.0.9.9000 (2015-03-25)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Changed the behaviour of \code{setConstraints} in that it now
  returns merely a warning instaed of an erro if more or less than the
  necessary constraints are specified.
    }
  }
}



\section{Changes in brRasch version 0.0.8.9000 (2015-03-25)}{
  \subsection{ENHANCEMENTS}{
    \itemize{
      \item Added support for \CRANpkg{Matrix package}
      \item Updated \code{?brRasch}
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Corrected a bug in \code{setConstraints} that was causing wrong
      constraints if the vector supplied in the \code{which} argument was
      not ordered.
    }
  }

}

\section{Changes in brRasch version 0.0.7.9000 (2015-03-24)}{

  \subsection{ENHANCEMENTS}{
    \itemize{
      \item Added support of \code{plot} method for \code{brRasch} objects with \code{dim = 2}
    }
  }

}


