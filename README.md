# OncoSignature - Drug response prediction on Acute Myeloid Leukemia

Copyright (C) 2019 Nicolàs Palacio<br>
Contact: nicolas.palacio@bioquant.uni-heidelberg.de<br>

GNU-GLPv3:

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

A full copy of the GNU General Public License can be found on
http://www.gnu.org/licenses/.

This repository contains the scripts used to analyze and extract
phosphorylation signatures of drug response in AML samples.

![Grapical abstract](oncosignature_pipeline.svg)

## 1) Data normalization pipeline
Data cleaning + normalization + batch correction

## 2) Differential expression analysis
This script analyzes the differential expression between treated vs. untreted
samples (responders and non-responders seperatedly) and between responders
vs. non-responders before treatment

### 2.1) Differential expression plots
This script generates the volcano plots to visualize the differential
expression results

### 2.2) Biomarkers of interest
This script tries to find interesting biomarkers based on the differential
expression results. Phosphosites of interest are defined as follows:

### 2.3) Significant differential phosphorylation
This script subsets the differential expresssion results and saves
a table for each contrast containing only the significantly
differentially expressed p-sites.

## 3) Gene Set Enrichment Analysis
This script computes the GSEA with 11 different methods and generates a
consensus score based on ranking

### 3.1) Gene set enrichment analysis plots
This script generates the plots to visualize the results of the GSEA

## 4) Kinase-substrate enrichment analysis
This script computes the kinase enrichment based on the differential
phosphoproteomics data
