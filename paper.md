**Transcriptionary: a Python library for interactive and customizable gene transcript visualizations**

Suchita Lulla<sup>1</sup>, Aaron Quinlan<sup>1</sup><br>
<sup>1</sup>Department of Human Genetics, University of Utah, Salt Lake City, UT

**Summary**<br>
Transcriptionary is a Python 3 command-line tool that accepts user-defined coordinate, variant, and annotation data in various file formats to create interactive and customizable visualizations on transcript models. It can integrate multiple data sets and annotations to visualize patterns of genetic variation with quantitative measures and annotations such as sequencing coverage and protein domain annotations. The visualizations from Transcriptionary can be saved as interactive, serverless HTML documents or exported to vector graphics formats.

**Statement of Need**<br>
Although various genome browsers exist, such as UCSC Genome Browser (Kuhn, Haussler and Kent 2013), IGV (Thorvaldsdóttir, Robinson and Mesirov 2013), and JBrowse (Buels et al. 2016), these tools depict entire transcripts and lack the functionality to remove or compress introns to maximize space for visualizing patterns on exons. Resources such as the gnomAD (Karczewski et al. 2021) browser provide this functionality but lack the ability to display user-defined data or customize the visualization. Moreover, MutationMapper (Jay and Brouwer 2016) uses cBioPortal (Gao et al. 2016) to create "lollipop" plots for mutations, but these plots are only created in protein space, not gene space. Transcriptionary provides flexibility in plotting user-supplied data along gene transcripts. It allows the user to control the screen space dedicated to intronic and exonic regions, and visualizations featuring user data and “lollipop” plots can be produced in PNG or SVG format or as a single, interactive HTML file.

**Usage**<br>
Transcriptionary uses the Python data visualization library bokeh (Bokeh Development Team 2019) to create gene transcript visualizations. To minimize the visual space attributed to introns, the library converts all user-supplied genomic coordinates to adjusted coordinates by compressing the introns to be 10 base pairs in length by default (Fig. 1a, 1b). User-supplied annotations, such as genetic variants and protein domains, are also converted to exonic coordinates (Fig. 1c). All plot elements are annotated with both genomic and exonic coordinates in an interactive hover box.

![alt text](https://github.com/quinlan-lab/transcriptionary/blob/main/transcriptionary_figure_final.jpg?raw=true)
**Figure 1. Overview of Transcriptionary's functionality:** (a) When visualizing an uncompressed transcript, introns visually take up almost all of the plot. (b) Introns are compressed to 10 bp to view the transcript more concisely. (c) Information from other datasets is likewise compressed and displayed alongside the exons, such as sequencing coverage shown by the gray bars, measures of genetic variation density shown by the blue and red lines, and protein domains shown by the green track.

Transcriptionary provides the ability to view multiple transcripts of the same gene on a single visualization, as some analyses necessitate plotting multiple isoforms next to each other. The user can specify multiple transcripts to plot via their transcript IDs, or they can choose to plot all transcripts for a given gene found in a gene annotation file in either GFF (gff3.md at Master · The-Sequence-Ontology/Specifications) or BED (Jeffrey Niu, Danielle Denisko, Michael M. Hoffman 2022) format. Additionally, the user can plot a "flattened" transcript, which combines all the transcripts for the gene found in the gene annotation file into a union of all exons. A flattened transcript enables, for example, visualizing the set of all variants that could be exonic in any isoform. 

Transcriptionary allows the user to visualize other annotations alongside one or more isoforms. The user can supply a VCF (Danecek et al. 2011) or BED file to display genetic variants, which are plotted as lollipops and automatically annotated with user-specified fields such as nucleotide change and variant allele frequency. If variant effect annotations from VEP (McLaren et al. 2016) are provided, variants are colored by severity and annotated with specified VEP fields. The user can choose to either use VEP annotations specific to each transcript or default to the annotation with the highest impact severity. Lollipops can be hidden and unhidden according to severity with a checkbox. User-provided track annotations, such as protein domains and genetic repeats, can be specified as GFF, GTF, or BED and are annotated with user-specified fields. The user can also supply coordinate data, such as sequencing coverage, as either TSV, CSV, or BEDGRAPH. The user can choose to include a slider that can smooth the coordinate data.
  
The y-axis for lollipop heights can be toggled using radio buttons to represent allele count or allele frequency, and the scale can be toggled to be linear or log. Checkboxes allow users to turn plot elements such as UTRs, domains, and regional data on and off. The colors of all plot elements can be specified individually, or the user can choose pre-defined, color-blind friendly color palettes (Tol; Wong 2011). The user can mouse over plot elements such as lollipops to view annotated data and click an exon to expand it to the entire width of the plot.

Transcriptionary visualizations can be saved as standalone, interactive HTML documents or static PNG or SVG images. Collectively, the Transcriptionary's options provide researchers with a flexible tool to create custom visualizations of gene transcripts and genomic annotations with compressed introns, allowing the user to focus their analyses on coding sequence.

**Installation**<br>
The Python packages necessary for running transcriptionary can be installed via a requirements.txt file, and transcriptionary can be installed through PyPI. Detailed instructions are available at the following URL: https://github.com/quinlan-lab/transcriptionary#install. To run transcriptionary, the user must provide a configuration file, and an example configuration file describing the parameters is available in the README file at the GitHub repository.

**Acknowledgments**<br>
The authors would like to thank Brent Pedersen for contribution to packaging the software and Laurel Hiatt for feedback on the manuscript.

**Funding**<br>
This research was funded by a grant (R01HG012252) from the NIH to ARQ.





**References**<br>
Bokeh Development Team. Bokeh: Python Library for Interactive Visualization., 2019.

Buels R, Yao E, Diesh CM et al. JBrowse: a dynamic web platform for genome visualization and analysis. Genome Biol 2016;17:66.

Danecek P, Auton A, Abecasis G et al. The variant call format and VCFtools. Bioinformatics 2011;27:2156–8.

Gao J, Lindsay J, Watt S et al. Abstract 5277: The cBioPortal for cancer genomics and its application in precision oncology. Cancer Res 2016;76:5277–5277.

gff3.md at Master · The-Sequence-Ontology/Specifications. Github

Jay JJ, Brouwer C. Lollipops in the Clinic: Information Dense Mutation Plots for Precision Medicine. PLoS One 2016;11:e0160519.

Jeffrey Niu, Danielle Denisko, Michael M. Hoffman. The Browser Extensible Data (BED) Format., 2022.

Karczewski KJ, Francioli LC, Tiao G et al. Author Correction: The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 2021;590:E53.

Kuhn RM, Haussler D, Kent WJ. The UCSC genome browser and associated tools. Brief Bioinform 2013;14:144–61.

McLaren W, Gil L, Hunt SE et al. The Ensembl Variant Effect Predictor. Genome Biol 2016;17:122.

Thorvaldsdóttir H, Robinson JT, Mesirov JP. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration. Brief Bioinform 2013;14:178–92.<br>

Tol P. Paul Tol’s Notes.

Wong B. Color blindness. Nat Methods 2011;8:441.

