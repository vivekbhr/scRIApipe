#!/usr/bin/env Rscript

## Usage: cat velocity_report.R | R --vanilla --quiet --args
.libPaths(R.home("library"))

args = commandArgs(TRUE)
rmdTemplate <- args[1]

library(magrittr)
library(BUSpaRse)
#library(velocyto.R)
library(Seurat)
## data transformation
library(ggplot2)
library(tibble)
library(dplyr)
library(uwot)
library(SeuratWrappers)
library(RColorBrewer)
library(AnnotationDbi)
library(rmarkdown)

file.copy(rmdTemplate, to = 'velocity_report.Rmd')
render('velocity_report.Rmd', output_format = "html_document", clean = TRUE)
