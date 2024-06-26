from .process_gene_gff import gff_to_db, get_gene_feature, get_transcript_dict
from .get_coords import get_variants,get_line,get_track
from .colors import color_boxes
from .axes import add_user_axis,add_variant_axis,add_axis
from .glyphs import add_intron_glyph, add_exon_glyph, add_variant_glyph, add_UTR_glyph, add_track_glyph, add_multi_line_glyph
from .widget_callbacks import add_checkbox,add_variant_severity_checkbox,add_user_tracks_checkbox,add_user_lines_checkbox,add_smoothing_slider,add_legend,add_linear_log_scale,add_exon_zoom,add_variant_sets_checkbox
from . import project_coords
import numpy as np
import argparse
from bokeh.plotting import figure, output_file, save
from bokeh.models.tools import HoverTool
from bokeh.layouts import column, gridplot
from bokeh.models import ColumnDataSource, Range1d, HoverTool, LabelSet, Div
import subprocess
import yaml
from yaml.loader import SafeLoader

def plot_transcript(plot_params, variant_params, user_track_params, user_line_params, transcript_dict, glyph_dict, variant_axes, line_axes, user_tracks, user_track_glyphs, user_lines, user_line_glyphs, title=''):
    project_coords.adjust_coordinates(transcript_dict['exons'], intron_size=plot_params['intron_size'])
    plot = figure(title=title, width=1500, tools='xpan,xzoom_in,xzoom_out,reset,save', height=plot_params['plot_height'],min_border=0,
               x_range=Range1d(0, max([e['compact_end'] for e in transcript_dict['exons']]), bounds='auto'), y_range=Range1d(0,plot_params['plot_height']), background_fill_color='white')
    plot.toolbar.logo = None
    plot.grid.grid_line_color = None
    plot.yaxis[0].visible = False

    transcript_ID = transcript_dict['ID'].split(':')[-1]
    
    ### VARIANTS ###

    for variant_set in variant_params: # need to map_point all variant sets first because add_variant_glyph uses 'compact_pos' from all variant sets
        project_coords.map_point(variant_params[variant_set]['variant_ls'], transcript_dict['exons'])

    for variant_set in variant_params:

        ray_glyph,circle_glyph = add_variant_glyph(plot_params, variant_params, variant_set, transcript_ID, plot)

        if ray_glyph and circle_glyph:
            tooltips_variant_set = [('Position (compact)', '@x'), ('Position (chr)', '@pos'), ('Severity', '@sev')]
            if variant_params[variant_set]['format'] == 'vcf': tooltips_variant_set.extend([('REF', '@ref'), ('ALT', '@alt')]) #for vcf automatically annotate REF and ALT
            tooltips_variant_set.extend([(info_field, '@'+info_field) for info_field in variant_params[variant_set]['info_annotations']]) #add user defined INFO fields (from vcf) to hover box annotation
            tooltips_variant_set.extend([(vep_field, '@'+vep_field) for vep_field in variant_params[variant_set]['vep']['vep_fields']]) #add user defined VEP field to hover box annotation
            if variant_params[variant_set]['has_yaxis_info']: tooltips_variant_set.extend([('Allele count', '@allele_counts'), ('Allele frequency', '@allele_frequencies')])

            if len(variant_params) > 1: #if more than one variant dataset, add dataset name to hover box annotation
                tooltips_variant_set = [('Dataset', '@variant_set')] + tooltips_variant_set

            plot.add_tools(HoverTool(tooltips=tooltips_variant_set, renderers=[circle_glyph,ray_glyph], point_policy='follow_mouse', attachment='below', visible=False))

            glyph_dict['Variant'].extend([ray_glyph,circle_glyph])

        else:
            print('Variant set \'{}\' has no variants that fall in transcript {}.'.format(variant_set, transcript_dict['ID']))

    if plot_params['add_variant_axis']:
        def log10(f): return np.log10(f) if f > 0 else 0

        all_vars = sum([variant_params[var_set]['variant_ls'] for var_set in variant_params], []) #flatten list 
        all_vars_in_transcript = [v for v in all_vars if v['compact_pos'] >= 0] #get list of all variants in transcript to get max and min for variant axes
        allele_counts = [v['allele_count'] for v in all_vars_in_transcript]
        allele_frequencies = [v['allele_frequency'] for v in all_vars_in_transcript]

        variant_axes['count_linear'].append(add_variant_axis(plot_params, variant_params, plot, 'Allele count', allele_counts, visible=(plot_params['default_y_axis'] == 'AC' and plot_params['default_y_axis_scale'] == 'linear')))
        variant_axes['count_log'].append(add_variant_axis(plot_params, variant_params, plot, 'Log(Allele count)', [log10(c) for c in allele_counts], visible=(plot_params['default_y_axis'] == 'AC' and plot_params['default_y_axis_scale'] == 'log')))
        variant_axes['frequency_linear'].append(add_variant_axis(plot_params, variant_params, plot, 'Allele frequency', allele_frequencies, visible=(plot_params['default_y_axis'] == 'AF' and plot_params['default_y_axis_scale'] == 'linear')))
        variant_axes['frequency_log'].append(add_variant_axis(plot_params, variant_params, plot, 'Log(Allele frequency)', [log10(f) for f in allele_frequencies], visible=(plot_params['default_y_axis'] == 'AF' and plot_params['default_y_axis_scale'] == 'log'))) #don't include 0s since they will just be plotted as min

    ######

    tooltips_features = [('Type','@feat_type'), ('Start (compact)', '@adj_start'), ('End (compact)', '@adj_end'), 
                      ('Start (chr)', '@true_start'), ('End (chr)', '@true_end'), ('Length', '@true_len')]

    ### INTRONS ###
    introns = project_coords.get_introns_from_exons(transcript_dict['exons'])
    introns_glyph = add_intron_glyph(plot_params, plot, introns)
    plot.add_tools(HoverTool(tooltips=tooltips_features, renderers=[introns_glyph], point_policy='follow_mouse', attachment='below', visible=False))    
    
    ### EXONS ###
    exons_glyph,arrow_glyph = add_exon_glyph(plot_params, plot, transcript_dict['exons'], transcript_dict['direction'])
    plot.add_tools(HoverTool(tooltips=tooltips_features, renderers=[exons_glyph], point_policy='follow_mouse', attachment='below', visible=False))
    
    glyph_dict['Direction'].append(arrow_glyph)
    glyph_dict['exon'].append(exons_glyph)

    ### UTRs ###
    project_coords.map_box(transcript_dict['UTRs'], transcript_dict['exons'])
    UTR_glyph = add_UTR_glyph(plot_params, plot, transcript_dict['UTRs'])
    plot.add_tools(HoverTool(tooltips=tooltips_features, renderers=[UTR_glyph], point_policy='follow_mouse', attachment='below', visible=False))    
    glyph_dict['UTRs'].append(UTR_glyph)

    ### USER TRACKS ###
    h = plot_params['track_height']

    for idx,track_name in enumerate(user_tracks):
        project_coords.map_box(user_tracks[track_name], transcript_dict['exons'])
        tooltips_tracks = [(user_track_params[track_name]['color_by'],'@color_by'), ('Start (adjusted)', '@adj_start'), ('End (adjusted)', '@adj_end'), 
                    ('Start (true)', '@true_start'), ('End (true)', '@true_end'), ('Length', '@true_len')]
        tooltips_tracks.extend([(field, '@'+field) for field in user_track_params[track_name]['annotate_with']])

        y = ((h*1.5)*(len(user_tracks) - idx - 1)+(h))
        track_glyph = add_track_glyph(user_track_params, track_name, plot, user_tracks[track_name], h*0.9, y)
        plot.add_tools(HoverTool(tooltips=tooltips_tracks, renderers=[track_glyph], point_policy='follow_mouse', attachment='below', visible=False))
        user_track_glyphs[track_name].append(track_glyph)
        cs = ColumnDataSource(dict(x=[0], y=[(h*1.5)*(len(user_tracks) - idx - 1)], text=[list(user_tracks.keys())[idx]]))
        label = LabelSet(source=cs, x='x', y='y', text='text',text_font_size='{}px'.format(plot_params['track_height']), text_align='left')
        plot.add_layout(label)
        user_track_glyphs[track_name].append(label)
        
    ### USER LINES ###
    for idx,axis_name in enumerate(user_line_params):
        all_xs = []
        all_ys = []
        for line in user_line_params[axis_name]['lines']:
                xs_ls,ys_ls = project_coords.map_line(user_lines[axis_name][line], transcript_dict['exons'], user_line_params[axis_name]['lines'][line]['chrom'])
                all_xs.append(xs_ls)
                all_ys.append(ys_ls)
        try: y_max = max([i for s in [i for s in all_ys for i in s] for i in s]) #flatten 3D list to 1D list to take max
        except: continue

        for idx,line in enumerate(user_line_params[axis_name]['lines']):
            line_params=user_line_params[axis_name]['lines'][line]
            line_glyph = add_multi_line_glyph(plot_params, plot, all_xs[idx], all_ys[idx], max_y=y_max, y0=plot_params['y0'], fill_area=line_params['fill_area'], line_color=line_params['color'], line_alpha=line_params['alpha'])
            line_glyph.level = 'underlay'
            user_line_glyphs[axis_name].append(line_glyph)
            
        line_axes[axis_name].append(add_user_axis(plot, plot_params, user_line_params, axis_name, y_max, plot_params['y0'], plot_params['plot_height'], y_min=0, num_ticks=3, axis_position='right', visible=True))

    return plot,glyph_dict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help='must be .yaml; sample config file can be found at test/test.yaml')
    args = parser.parse_args()

    ### CONFIG ###
    config_file = args.config_file
    if config_file[-5:] != '.yaml': raise RuntimeError('configuration file must be .yaml')
    with open(config_file) as f:
        params = list(yaml.load_all(f, Loader=SafeLoader))
        plot_params = params[0]
        variant_params = params[1]
        user_track_params = params[2]
        user_line_params = params[3]
        if not variant_params: variant_params = []
        if not user_track_params: user_track_params = []
        if not user_line_params: user_line_params = []

        plot_params['transcript_height'] = 20 + len(user_track_params) * plot_params['track_height'] * 1.5
        plot_params['y0'] = plot_params['transcript_height'] + plot_params['exon_height'] / 2  # "0" for line plots

    ### COLORS ###
    with open(plot_params['named_colors_filepath']) as f:
        named_colors = list(yaml.load_all(f, Loader=SafeLoader))[0]

    ### GLYPH COLORS ###
    # replace named colors with hex code from named_colors.yaml
    def get_color(color): return color if color[0] == "#" else named_colors[color]

    for glyph_type in plot_params['glyph_colors']: plot_params['glyph_colors'][glyph_type] = get_color(plot_params['glyph_colors'][glyph_type])

    for variant_set in variant_params:
        for sev in variant_params[variant_set]['variant_severity_colors']:
            variant_params[variant_set]['variant_severity_colors'][sev] = get_color(variant_params[variant_set]['variant_severity_colors'][sev])
        variant_params[variant_set]['color'] = get_color(variant_params[variant_set]['color'])

    for track in user_track_params:
        try:
            for domain in user_track_params[track]['colors']:
                user_track_params[track]['colors'][domain] = get_color(user_track_params[track]['colors'][domain])
        except: continue

    with open(plot_params['palettes_filepath']) as f:
        palettes = list(yaml.load_all(f, Loader=SafeLoader))[0]
    track_colors = [c if '#' in c else named_colors[c] for c in palettes[plot_params['track_palette']]]
    plot_params['track_colors'] = track_colors

    for axis in user_line_params:
        for line in user_line_params[axis]['lines']:
            user_line_params[axis]['lines'][line]['color'] = get_color(user_line_params[axis]['lines'][line]['color'])
    
    ### VARIANTS ###
    for variant_set in variant_params:

        variant_params[variant_set]['format'] = variant_params[variant_set]['format'].lower().strip('.')

        try: #if info_annotations empty or nonexistent, set to empty list
            if not variant_params[variant_set]['info_annotations']: 
                variant_params[variant_set]['info_annotations'] = []
        except: variant_params[variant_set]['info_annotations'] = []

        try: #if consequence_idx nonexistent,
            int(variant_params[variant_set]['consequence_idx'])
        except: variant_params[variant_set]['consequence_idx'] = False

        #for BED files, replace space with underscore in annotation field names (space causes issues with hover box)
        try:
            variant_params[variant_set]['header'] = list(map(lambda field_name: str.replace(field_name, ' ', '_'), variant_params[variant_set]['header']))
            variant_params[variant_set]['info_annotations'] = list(map(lambda field_name: str.replace(field_name, ' ', '_'), variant_params[variant_set]['info_annotations']))
        except: pass

        try: #if vep empty or nonexistent, set to empty params
            if not variant_params[variant_set]['vep']: variant_params[variant_set]['vep'] = {'field_name': '' ,'vep_fields': [], 'annotate_severity_by': ''}
        except: variant_params[variant_set]['vep'] = {'field_name': '' ,'vep_fields': [], 'annotate_severity_by': ''}

    ### TRACKS ###
    for track_name in user_track_params:
        try: #if annotate_with empty or nonexistent, set to empty list
            if not user_track_params[track_name]['annotate_with']: user_track_params[track_name]['annotate_with'] = []
        except: user_track_params[track_name]['annotate_with'] = []

    ### LINES ###
    if not user_line_params: user_line_params = []

    ### TRANSCRIPTS ###
    transcript_IDs = plot_params['transcripts']

    ### OUTPUT ###
    output_format = plot_params['output_format'].lower().strip('.')
    if output_format not in ['html', 'png', 'svg']: raise RuntimeError('Argument output_format in {} must be HTML, PNG, or SVG'.format(config_file))
    output = plot_params['output_filepath'] + '.' + output_format if plot_params['output_filepath'][-len(output_format):] != output_format else plot_params['output_filepath']

    return plot_params,variant_params,user_track_params,user_line_params,transcript_IDs,output,output_format

def transcriptionary():
    plot_params, variant_params, user_track_params, user_line_params, transcript_IDs, output, output_format = parse_args()
    gff_db = gff_to_db(plot_params['gff_path'],plot_params['gff_path']+'.db')
    
    gene_feature = get_gene_feature(gff_db, plot_params['gene_name'])

    if transcript_IDs == 'transcript_names':
        transcripts = get_transcript_dict(plot_params, gff_db, gene_feature, 'all')
        print(list(transcripts.keys()))
        exit()

    transcripts = get_transcript_dict(plot_params, gff_db, gene_feature, transcript_IDs)
    transcript_IDs = list(transcripts.keys()) #if 'all', transcript_IDs will become list of transcript names; if nonexistent IDs they are removed

    plot_params['plot_variants'] = bool(variant_params)
    plot_params['add_variant_severity_checkbox'] = False #will be set to true if consequence annotations are present

    for variant_set in variant_params:
        variant_params[variant_set]['variant_ls'] = get_variants(plot_params, variant_params, variant_set, transcripts, gene_feature.start, gene_feature.end) if plot_params['plot_variants'] else []

    plot_params['add_variant_axis'] = len([variant_set for variant_set in variant_params if variant_params[variant_set]['has_yaxis_info']]) >= 1

    user_tracks = {track_name: get_track(user_track_params, track_name) for track_name in user_track_params}
    for track_name in user_tracks: color_boxes(plot_params, user_track_params, track_name, user_tracks[track_name])
    user_lines = {axis_name:{} for axis_name in user_line_params}
    for axis_name in user_line_params:
        for line_name in user_line_params[axis_name]['lines']:
            user_lines[axis_name][line_name] = get_line(user_line_params, axis_name, line_name)
    
    plot_ls = []
    glyph_dict = dict(exon=[],UTRs=[],Variant=[],Direction=[])
    user_track_glyphs = {track_name:[] for track_name in user_track_params}
    user_line_glyphs = {line_name:[] for line_name in user_line_params}
    variant_axes = dict(count_linear=[],count_log=[],frequency_linear=[],frequency_log=[])
    line_axes = {}
    for line_name in user_line_params:
        line_axes[line_name] = []
    
    for ID in transcript_IDs:
        title = 'gene={}; transcript={}/{}'.format(plot_params['gene_name'], ID, transcripts[ID]['ID'])
        if transcripts[ID]['direction']: title += ' ({})'.format(transcripts[ID]['direction'])
        plot,glyph_dict = plot_transcript(plot_params, variant_params, user_track_params, user_line_params, transcripts[ID], glyph_dict, variant_axes, line_axes, user_tracks, user_track_glyphs, user_lines, user_line_glyphs, title=title)
        plot_ls.append(plot)

    legend = [add_legend(user_line_params)] if user_line_params else []

    if output_format == 'html': #only add widgets for HTML  

        empty_plot = figure(height=1500,outline_line_color=None, toolbar_location=None) #for HTML, add white space at bottom so hover boxes are not cut off 
        empty_plot.line(x=[0], y=[0]) #avoid empty plot warning
        empty_plot.yaxis.visible = empty_plot.xaxis.visible = empty_plot.grid.visible = False
        # plot_ls.append(empty_plot) 

        add_exon_zoom(plot_ls,glyph_dict)
        checkbox = add_checkbox(glyph_dict,plot_params)
        
        user_line_checkboxes=[] #one checkbox per user line, so that they can be lined up with sliders
        for axis in user_line_params:
            user_line_checkbox = add_user_lines_checkbox(plot_ls, line_axes[axis], user_line_glyphs[axis], axis)
            user_line_checkboxes.append(user_line_checkbox)
        
        sliders = []
        for axis_name in user_line_params:
            if user_line_params[axis_name]['smoothing_slider']:
                fill_area_ls = [user_line_params[axis_name]['lines'][line]['fill_area'] for line in user_line_params[axis_name]['lines']]*len(plot_ls)
                sliders.append(add_smoothing_slider(user_line_glyphs[axis_name], fill_area_ls, title='{} smoothing'.format(axis_name)))
            else:
                sliders.append(None)

        grid1 = [[checkbox]]
        if plot_params['plot_variants']:
            if plot_params['add_variant_severity_checkbox']:
                variant_severity_checkbox = add_variant_severity_checkbox(plot_ls, glyph_dict)
                grid1[0].append(variant_severity_checkbox)
            if len(variant_params) > 1: #if more than one variant_set to plot
                variant_sets_checkbox = add_variant_sets_checkbox(glyph_dict, variant_params)
                grid1[0].append(variant_sets_checkbox)
            if plot_params['add_variant_axis']:
                extra_empty_divs_1,extra_empty_divs_2 = ([Div(text=""" """, width=200, height=15) for i in range(len(grid1[0]) - 2)] for j in range(2))
                div_type,radio_group_type,div_scale,radio_group_scale = add_linear_log_scale(plot_params, variant_axes, glyph_dict)
                grid1.extend([[div_type, div_scale] + extra_empty_divs_1, [radio_group_type, radio_group_scale] + extra_empty_divs_2])

        if len(user_track_glyphs) > 0:
            axes = {**variant_axes,**line_axes}

            if not [ax_ls for ax_ls in axes.values() if ax_ls]: #if no user-defined y axis, add empty invisible y axis so tracks can move according to checkbox
                axes = {'empty': [add_axis(plot, plot_params, 'empty_axis', 1, 0, plot_params['plot_height'], 1, False, visible=False)]}

            user_tracks_checkbox = add_user_tracks_checkbox(plot_ls,axes,user_track_glyphs,glyph_dict['Direction'],plot_params)
            grid1[0].append(user_tracks_checkbox)

        lines = list(zip(user_line_checkboxes,sliders))
        for tup in lines: grid1.append(tup)
        grid1.append(legend)
        grid = gridplot(grid1, merge_tools=False)
        output_file(output)
        save(column([grid]+plot_ls))

        #change HTML title from 'Bokeh Plot'
        subprocess.run(['sed', '-i', 's/<title>Bokeh Plot<\/title>/<title>{}<\/title>/g'.format(plot_params['title']), output])

    elif output_format == 'png':
        from bokeh.io import export_png
        for p in plot_ls: p.toolbar_location = None
        export_png(column(legend + plot_ls), filename=output)

    elif output_format == 'svg':
        from bokeh.io import export_svg
        for p in plot_ls: p.output_backend = 'svg'
        export_svg(column(legend + plot_ls), filename=output)

main = transcriptionary
