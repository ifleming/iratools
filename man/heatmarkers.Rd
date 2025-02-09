% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmarkers.R
\name{heatmarkers}
\alias{heatmarkers}
\title{Generate a heatmap of scRNA-seq differential expression results}
\usage{
heatmarkers(
  seurat_object,
  markers_df,
  n_markers_shown,
  identity_comparison,
  slot,
  ranking_metric = "avg_log2FC",
  cluster = F
)
}
\arguments{
\item{seurat_object}{Seurat object containing the single-cell data}

\item{markers_df}{Output dataframe from Seurat's FindAllMarkers() containing genes, logfold change, and p-values}

\item{n_markers_shown}{Number of positive markers per condition to show in heatmap}

\item{identity_comparison}{Column name in Seurat object metadata used to group cells e.g. "cell_type"}

\item{slot}{Which data slot to map}

\item{ranking_metric}{Column name to rank genes by (default: "avg_log2FC")}

\item{cluster}{Decision on whether genes should be clustered or just listed by gene rank}
}
\value{
A ggplot object containing the heatmap
}
\description{
Generate a heatmap of scRNA-seq differential expression results
}
\author{
Ira Fleming, \email{ira.fleming@cuanschutz.edu}
}
