# transcriptionary: customizable, interactive gene transcript plots

transcriptionary takes user-defined parameters to create a static .html document displaying gene transcripts with introns compressed. 

Exons are annotated with coordinates and length.

Variants are specified as .vcf, .csv, or .bed and are annotated in hover boxes with coordinate, allele count, allele frequency, mutation type, and VEP annotation. They are colored by severity. Radio buttons in the plot give the option to view lollipop heights as allele count or allele frequency, and linear or log scale. 
If a VCF with VEP annotations is provided, lollipops will be categorized as LOW, MED, or HIGH impact as per [geneimpacts](https://github.com/brentp/geneimpacts) and colored accordingly.

Tracks are specified as .gtf and are annotated with name, coordinates, and length. Colors can be specified by the user or chosen randomly from a color palette.

Coordinate-based information can be provided as .csv/.tsv (point-based) or .bedgraph (interval-based). The user can customize the y axis with tick precision and scientific notation. The user can specify the line color, alpha value, and choose whether to fill in the area under the curve.

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

`gff_path`: path to gff (for feature coordinates). When running the first time, a `gff.db` file will be created for you. When rerunning, can change this parameter to the .gff.db file to avoid recreating it.

`variant_path`: path to vcf, csv, or bed (for variant coordinates)

`gene_name`: gene name

`seqid`: chromosome

`plot_height`: height of plot in pixels; default 200
`plot_width`: width of plot in pixels; default 1500
`track_height`: height of tracks in pixels; default 10
`exon_height`: height of exons in pixels; default 16
`intron_size`: size to which introns are compressed; default 20

`plot_UTRs`: show UTRs (boolean)
`plot_direction`: show arrows with direction (boolean)

`glyph_colors`: specify feature colors (can be hex code or any name from default_colors/named_colors.yaml)
- `intron`: default 'gray 11'
- `exon`: default 'davys gray'
- `arrow`: default '#252525'
- `UTR`: default '#969696'
- `variant`: default lollipop color if no severity information is available; default 'charcoal'

`track_palette`: palette to draw random track colors from; can be any palette in default_colors/palettes.yaml

`plot_variants`: show lollipops (boolean)
`seqid`: chromosome
`min_lollipop_height`: minimum height of lollipop in pixels; default 15
`lollipop_radius`: radius of lollipop in pixels; default 5
`lollipop_line_width`: line width of lollipop in pixels; default 2
`variant_severity_colors`: specify lollipop colors by variant severity; use hex codes or predefined colors from default_colors/named_colors.yaml
- LOW: default '#80b918'
- MODERATE: default '#f6aa1c'
- HIGH: default '#e01e37'
- MODIFIER: default '#778da9'

`track_palette`: any palette in default_colors/palettes.yaml

`<track_name>`: track label
- `gtf_path`: path to .gtf with track coordinate information 
- `seqid`: chromosome
- `colors`: specify colors by domain name; if no color is specified for a domain name, it will be chosen from the palette specified in track_palette
      -`<domain_name>`: (can be hex code or any name from default_colors/named_colors.yaml)

`<axis_name>`:
- `y_axis_label`: y axis label
- `num_ticks`: number of y ticks; default 3
- `tick_precision`: number of decimal places
- `tick_scientific_notation`: (bool)
- `smoothing slider`: include slider widget to smooth lines (boolean)
- `lines`: lines to plot on this axis
    - `<line_name>`:
        - `filepath`: (.csv or .bedgraph)
        - `color`: 
        - `alpha`:
        - `fill_area`: fill area under the curve (boolean)


## Run Test
```
python scripts/transcriptionary.py test/test.yaml all --output test/plot.html
```

## How to customize colors and color palettes
Define or adjust named colors by modifying `default_colors/named_colors.yaml`. Create a custom color palette to pull random track colors from by adding to `default_colors/palettes.yaml` (use hex codes or use the named colors defined in `default_colors/named_colors.yaml`).

## Color palette citations
Color-blind friendly color palettes are included in default_colors/palettes.yaml

- bang_wong_palette: https://www.nature.com/articles/nmeth.1618
- paul_tol_palette: https://personal.sron.nl/~pault/