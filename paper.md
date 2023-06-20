---
title: 'Transcriptionary: a Python library for interactive and customizable gene transcript visualizations'
tags:
  - Python
  - gene
  - transcript
  - visualization
authors:
  - name: Suchita Lulla
    orcid: 0009-0007-2301-2951
    equal-contrib: true
    affiliation: 1
  - name: Aaron Quinlan
    orcid: 0000-0003-1756-0859
    equal-contrib: false
    affiliation: 1
    corresponding: true
affiliations:
  - name: University of Utah, USA
    index: 1
date: 20 June 2023
bibliography: paper.bib

---

# Summary

`Transcriptionary` is a Python 3 command-line tool that accepts user-defined coordinate, variant, and annotation data in various file formats to create interactive and customizable visualizations on transcript models. It can integrate multiple data sets and annotations to visualize patterns of genetic variation with quantitative measures and annotations such as sequencing coverage and protein domain annotations. The visualizations from `Transcriptionary` can be saved as interactive, serverless HTML documents or exported to vector graphics formats.

# Statement of need

Although various genome browsers exist, such as UCSC Genome Browser [@ucscgenome], IGV [@igv], and JBrowse [@jbrowse], these tools depict entire genes and lack the functionality to remove or compress introns to maximize space for visualizing patterns on exons. Resources such as the gnomAD [@gnomad] browser provide this functionality but lack the ability to visualize custom annotations and other user-defined data such as sequencing coverage along with mutations. This limitation also exists in cBioPortal’s MutationMapper [@cbioportal] and Lollipops [@lollipops], which create "lollipop" plots for mutations only using protein coordinates, not gene coordinates ("gene space"). Therefore, there is a need for an interactive tool to explore mutations and other relevant data in gene space to facilitate variant analysis and interpretation. Transcriptionary provides flexibility in plotting user-supplied data along gene transcripts, while allowing the user to adjust introns sizes to focus limited visualization space on exons. It allows the user to control the screen space dedicated to intronic and exonic regions, and visualizations featuring user data and “lollipop” plots can be produced in PNG or SVG format or as a single, interactive HTML file.

# Usage

`Transcriptionary` uses the Python data visualization library `bokeh` [@bokeh] to create gene transcript visualizations. To minimize the visual space attributed to introns, the library converts all user-supplied genomic coordinates to adjusted coordinates by compressing the introns to be 10 base pairs in length by default (Fig. 1a, 1b). User-supplied annotations, such as genetic variants and protein domains, are also converted to exonic coordinates (Fig. 1c). All plot elements are annotated with both genomic and exonic coordinates in an interactive hover box.

![Overview of `Transcriptionary`'s functionality: (a) When visualizing an uncompressed transcript, introns visually take up almost all of the plot. (b) Introns are compressed to 10 bp to view the transcript more concisely. (c) Information from other datasets is likewise compressed and displayed alongside the exons, such as sequencing coverage shown by the gray bars, measures of genetic variation density shown by the blue and red lines, and protein domains shown by the green track.\label{fig:overview}](https://github.com/quinlan-lab/transcriptionary/blob/main/transcriptionary_figure_final.jpg?raw=true)

`Transcriptionary` provides the ability to view multiple transcripts of the same gene on a single visualization, as some analyses necessitate plotting multiple isoforms next to each other. The user can specify multiple transcripts to plot via their transcript IDs, or they can choose to plot all transcripts for a given gene found in a gene annotation file in either GFF [@gff] or BED [@bed] format. Additionally, the user can plot a "flattened" transcript, which combines all the transcripts for the gene found in the gene annotation file into a union of all exons. A flattened transcript enables, for example, visualizing the set of all variants that could be exonic in any isoform. 

`Transcriptionary` allows the user to visualize other annotations alongside one or more isoforms. The user can supply a VCF [@vcf] or BED file to display genetic variants, which are plotted as lollipops and automatically annotated with user-specified fields such as nucleotide change and variant allele frequency. If variant effect annotations from VEP [@vep] are provided, variants are colored by severity and annotated with specified VEP fields. The user can choose to either use VEP annotations specific to each transcript or default to the annotation with the highest impact severity. Lollipops can be hidden and unhidden according to severity with a checkbox. User-provided track annotations, such as protein domains and genetic repeats, can be specified as GFF, GTF, or BED and are annotated with user-specified fields. The user can also supply coordinate data, such as sequencing coverage, as either TSV, CSV, or BEDGRAPH. The user can choose to include a slider that can smooth the coordinate data.
  
The y-axis for lollipop heights can be toggled using radio buttons to represent allele count or allele frequency, and the scale can be toggled to be linear or log. Checkboxes allow users to turn plot elements such as UTRs, domains, and regional data on and off. The colors of all plot elements can be specified individually, or the user can choose pre-defined, color-blind friendly color palettes [@tolcb][@wongcb]. The user can mouse over plot elements such as lollipops to view annotated data and click an exon to expand it to the entire width of the plot.

`Transcriptionary` visualizations can be saved as standalone, interactive HTML documents or static PNG or SVG images. Collectively, the Transcriptionary's options provide researchers with a flexible tool to create custom visualizations of gene transcripts and genomic annotations with compressed introns, allowing the user to focus their analyses on coding sequence.

# Installation

The Python packages necessary for running `Transcriptionary` can be installed via a requirements.txt file, and `Transcriptionary` can be installed through PyPI. Detailed instructions are available at the following URL: https://github.com/quinlan-lab/transcriptionary#install. To run `Transcriptionary`, the user must provide a configuration file, and an example configuration file describing the parameters is available in the README file at the GitHub repository.

# Acknowledgments
The authors would like to thank Brent Pedersen for contribution to packaging the software and Laurel Hiatt for feedback on the manuscript.

# Funding
This research was funded by a grant (R01HG012252) from the NIH to ARQ.

# References
