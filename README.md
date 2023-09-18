[![workflow badge](https://github.com/quinlan-lab/transcriptionary/actions/workflows/python-package.yml/badge.svg)](https://github.com/quinlan-lab/transcriptionary/actions/workflows/python-package.yml)

# transcriptionary: customizable, interactive gene transcript plots

`transcriptionary`` takes user-defined parameters to create a static .html document displaying gene transcripts with introns compressed. 

Exons are annotated with coordinates and length.

Variants are specified as VCF or BED and are annotated with interactive hover boxes. For VCF, lollipops are automatically annotated with REF, ALT, allele count, allele frequency and are colored by severity if this information is provided. Radio buttons in the HTML plot give the option to view lollipop heights as allele count or allele frequency, and linear or log scale. 
If a variant file with VEP annotations is provided, lollipops will be categorized as LOW, MED, or HIGH impact as per [geneimpacts](https://github.com/brentp/geneimpacts) and colored accordingly. If no VEP annotation is present for a given transcript, lollipops on that transcript will be annotated as NONE. For a variant file with VEP annotations, lollipops can be turned on and off by severity with a checkbox. If more than one variant file is provided, each lollipop is annotated with the variant set it comes from, and each set of variants can be turned on and off with a checkbox.

Tracks are specified as GTF or BED and are annotated with name, coordinates, and length, along with other specified fields from the GTF/BED. Colors can be specified by the user or chosen randomly from a color palette.

Coordinate-based information can be provided as CSV/TSV (point-based) or BEDGRAPH (interval-based). The user can customize the y axis with tick precision and scientific notation. The user can specify the line color, alpha value, and choose whether to fill in the area under the curve.

If BED file contains a header, the line must begin with `#`. BED files must have at least three columns, which should be `chrom`, `start`, and `end`.

Plots can be output as HTML, PNG, or SVG.

For HTML plots: UTRs, direction arrows, and lollipops, tracks, and coordinate-based information can be turned on and off in the plot with checkboxes. A smoothing slider can be added to coordinate-based information. Click on any exon to expand it. Click on any white space to revert the plot.

## Demo Plot
https://home.chpc.utah.edu/~u6038618/transcriptionary/plot.html

## Arguments

`config_file` (required): user-defined parameters (sample config file at test/test.yaml)

## Config Parameters

`output_format`: can be HTML (interactive); PNG or SVG (non-interactive)

`output_filepath`: path to desired output file

`transcripts`: list transcript names or specify transcripts to plot. can be:
- `transcript-names`: list all possible transcript names for the given configuration file (does not create a plot).
- `flattened-exons`: overlay all transcripts to create a transcript with the largest possible exons (used to view all possible exonic variants).
- `all`: plot all transcripts, including `flattened-exons`.
- '[`<transcript_name_1>`, `<transcript_name_2>`, ...]': specify transcripts

`gff_path`: path to GFF or GTF (for feature coordinates). When running the first time, a `gff.db` file will be created for you. When rerunning, can change this parameter to the .gff.db file to avoid recreating it.

`gene_name`: gene name

`chrom`: chromosome (from features GFF/GTF)

`plot_height`: height of plot in pixels; default 200
`plot_width`: width of plot in pixels; default 1500
`track_height`: height of tracks in pixels; default 10
`exon_height`: height of exons in pixels; default 16
`intron_size`: size to which introns are compressed; default 10

`plot_UTRs`: show UTRs (boolean)
`plot_direction`: show arrows with direction (boolean)
`plot_variants`: show lollipops (boolean)

`min_lollipop_height`: minimum height of lollipop in pixels; default 15
`lollipop_radius`: radius of lollipop in pixels; default 5
`lollipop_line_width`: line width of lollipop in pixels; default 2

`default_y_axis`: set lollipop heights according to allele count (`AC`) or allele frequency (`AF`) by default (can be toggled with HTML output)
`default_y_axis_scale`: scale lollipop heights according on a linear (`linear`) or log (`log`) scale by default (can be toggled with HTML output)

`glyph_colors`: specify feature colors (can be hex code or any name from default_colors/named_colors.yaml)
- `intron`: default 'gray 11'
- `exon`: default 'davys gray'
- `arrow`: default '#252525'
- `UTR`: default '#969696'

`palettes_filepath`: filepath to config file with color palettes; default `default_colors/palettes.yaml`
`named_colors_filepath`: filepath to config file mapping hex color codes to named colors; default `default_colors/named_colors.yaml`
`track_palette`: palette to draw random track colors from; can be any palette in default_colors/palettes.yaml

`#VCF`
`<variant_set>`: variant set label
- `format`: file format of variant file; 'vcf'
- `filepath`: path to VCF
- `chrom`: chromosome
- `info_annotations`: for VCF, INFO fields to add to hover boxes
    - `<info_field_1>`
- `vep`: VCF only; leave empty if not VEP annotated
    - `field_name`: name of INFO field with VEP string (e.g. vep, ann, csq)
    - `vep_fields`: vep fields to add to hover boxes
        - `<vep_field_1>`
- `annotate_severity_by`: VCF only; possible arguments are `transcript_severity` (use VEP annotations from given transcript) and `max_severity` (use VEP annotation from transcript with most severe consequenced specified in `transcripts` argument); these apply to `vep_fields` also
- `color`: default lollipop color; use hex codes or predefined colors from default_colors/named_colors.yaml
- `variant_severity_colors`: specify lollipop colors by variant severity; use hex codes or predefined colors from default_colors/named_colors.yaml
    - `LOW`:
    - `MED`:
    - `HIGH`:

`#BED (if header in file, it must start with #)`
`<variant_set>`: variant set label
- `format`: file format of variant file; 'vcf'
- `filepath`: path to VCF
- `chrom`: chromosome
- `header`: BED only; list of column names in BED file
- `info_annotations`: for BED, column names to add to hover boxes
        - `<info_field_1>`
- `color`: default lollipop color; use hex codes or predefined colors from default_colors/named_colors.yaml
- `variant_severity_colors`: specify lollipop colors by variant severity; use hex codes or predefined colors from default_colors/named_colors.yaml
    - `LOW`:
    - `MED`:
    - `HIGH`:
- `consequence_idx`: BED only; index of file containing 'Consequence' field from VEP (leave empty or set to False if no such column); 0-based

`<track_name>`: track label
- `format`: file format of track file, can be GTF or BED
- `filepath`: path to file with track coordinate information 
- `seqid`: chromosome
- `header`: list of column names (BED only)
- `color_by`: field to color boxes by
- `annotate_with`: fields to annotate with
    - `<field_name>`: (for BED, must be from header)

`<axis_name>`:
- `y_axis_label`: y axis label
- `num_ticks`: number of y ticks; default 3
- `tick_precision`: number of decimal places
- `tick_scientific_notation`: (boolean)
- `smoothing slider`: include slider widget to smooth lines (boolean)
- `lines`: lines to plot on this axis
    - `<line_name>`:
        - `filepath`: (CSV or BEDGRAPH)
        - `color`: 
        - `alpha`:
        - `fill_area`: fill area under the curve (boolean)

## Install
```
git clone https://github.com/quinlan-lab/transcriptionary.git
cd transcriptionary
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n transcriptionary --file requirements.txt python=3.10
conda activate transcriptionary
python setup.py install
```

or

```
git clone https://github.com/quinlan-lab/transcriptionary.git
cd transcriptionary
pip install -r requirements.txt
python setup.py install
```

## Run Test
```
transcriptionary test/test.yaml
```

will create `transcriptionary-example.html`. Open `transcriptionary-example.html` in browser to view plots.

## How to customize colors and color palettes
Define or adjust named colors by modifying `default_colors/named_colors.yaml`. Create a custom color palette to pull random track colors from by adding to `default_colors/palettes.yaml` (use hex codes or use the named colors defined in `default_colors/named_colors.yaml`). Alternatively, create your own named_colors and palettes config files and change the file path in the main config file.

## Color palette citations
Color-blind friendly color palettes are included in default_colors/palettes.yaml

- bang_wong_palette: https://www.nature.com/articles/nmeth.1618
- paul_tol_palette: https://personal.sron.nl/~pault/
