% -*- coding: utf-8 -*-
\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options: width=0.5in}}}
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\sspace}{\ifelse{latex}{\out{~}}{ }}
%% some pre-defined commands: \R, \code, \acronym, \url, \file, \pkg

\name{TODO}
\title{Todo for package \CRANpkg{brRasch}}
\encoding{UTF-8}

\section{TODO after brRasch version 0.1.0.9000 (2015-03-29)}{
  \subsection{HIGH PRIORITY}{
    \itemize{
      \item Features: Add support for Jeffreys prior penalization
      \item Features: Add L1 and SCAD penalization capabilities
      \item Features: Add support for summary method
      \item Algorithm: Get the minimum of all steps if fails
      \item Algorithm: Handle NA likelihood values due to 0/1
      probabilities
      \item Algorithm: Check whether there can be an implicit stepsize
      calculation at the starting value level
      \item Algorithm: if very close to solution and ml then reduce the
      fsridge size
      \item Algorithm: control if step size is massive
      \item Algorithm: how to decide what Matrix class is used
      \item Use sparse or dense matrices?
    }
  }
  \subsection{LOW PRIORITY}{
    \itemize{
      \item Features: Add facilities to reshape from a subject-item  arrangement to a data frame ready for glm or gnm and back
      \item Features: Add facilities to move between scaled constrasts
      \item Features: Add item-response surfaces
      \item Features: Add chi-squared/likelihood ratio test facilities
    }
  }
}





