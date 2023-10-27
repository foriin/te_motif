# te_motif

[Shiny app](https://thenotoriousmbg.shinyapps.io/te_motif/) for screening for core promoter occurrences close to expressed transposons TSSs.


## Overview
Based on the CAGE data, I identified transcription start sites of active transposable elements (TEs) in Drosophila ovaries (ovaries being the main battlefield of genomic conflict between fly genome and transposons) in 3 conditions.Each condition reflects knockdown of one of the piRNA-pathway key protein. Then, based on the PFMs of Drosophila core promoter elements (CPEs) that were retrieved from [ElemeNT](https://www.juven-gershonlab.org/resources/element/) software (CLI application was kindly provided by Prof. Tamar Juven-Gershon upon request) I analyzed the occurences of CPEs close to TE TSSs (+- 40 nt). Stacked barplots show motif occurences for each nucleotide position.

## Usage
* Set the TPM range to select for TEs with high or low expression
* In TE selection tab you can select all, deselect all, or select specific Drosophila TE classes (LTR, LINE and DNA)
* In Motif selection tab you can include or remove additional CPEs
