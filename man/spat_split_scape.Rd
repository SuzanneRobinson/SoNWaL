% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialFuncsProjections.R
\name{spat_split_scape}
\alias{spat_split_scape}
\title{spat_split_scape}
\usage{
spat_split_scape(
  data_dir,
  save_file,
  output_dir,
  dates = c(1960:1990),
  variable = "all",
  num_chunks = c(1:10)
)
}
\arguments{
\item{data_dir}{directory which stores the CHESS spatial
data (under a filename called CHESS) - also should include
a folder for the full time series to be written into called fullTS}

\item{save_file}{where to save the output files}

\item{variable}{whether to process all CHESS climate variables or select
individual ones to process (e.g. precip, temp etc. )}

\item{num_chunks}{function creates chunks of 10,000
grid cells going down the UK, 1:39 covers scotland,
higher values cover the rest of UK}

\item{startDate}{the year to start the processing from, CHESS goes
back to 1961 but you may not need that far back
(significantly increases time to run function the further back you go)}
}
\value{
climate data with coordinates
}
\description{
spatial splitting function - split spatial
CHESS-scape data into chunks
}
