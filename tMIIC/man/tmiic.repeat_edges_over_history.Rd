% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmiic.wrapper.R
\name{tmiic.repeat_edges_over_history}
\alias{tmiic.repeat_edges_over_history}
\title{tmiic.repeat_edges_over_history}
\usage{
tmiic.repeat_edges_over_history(tmiic.res)
}
\arguments{
\item{tmiic.res}{[a tmiic object] The object returned by miic's 
execution in temporal mode.}
}
\value{
[a tmiic object] The tmiic object with a modified summary
}
\description{
Duplicates edges found by miic over the history assuming stationnarity
}
\details{
In temporal mode, the network returned by miic contains only edges
with at least one contemporaneous node (lag0). To improve the visual 
aspect when plotting, this function duplicates the edges over the history.  
i.e: assuming that we used tau=3, the edge X_lag0-X_lag1 will be 
copied as X_lag1-X_lag2 and X_lag2-X_lag3.\cr
Note that only the summary data frame is modified.
}
