from .process_gene_gff import gff_to_db, get_gene_feature, get_transcript_dict
from .get_coords import get_variants,get_line,get_track
from .colors import color_boxes
from .axes import add_user_axis,add_variant_axis
from .glyphs import add_intron_glyph, add_exon_glyph, add_variant_glyph, add_UTR_glyph, add_track_glyph, add_multi_line_glyph
from .widget_callbacks import add_checkbox,add_user_tracks_checkbox,add_user_lines_checkbox,add_smoothing_slider,add_legend,add_linear_log_scale,add_exon_zoom
from . import project_coords
import numpy as np
import argparse
from bokeh.plotting import figure#, output_file, save
from bokeh.layouts import column, row, gridplot
from bokeh.models import ColumnDataSource, Range1d, HoverTool, LabelSet
import yaml
from yaml.loader import SafeLoader

def plot_transcript(plot_params, variant_params, user_track_params, user_line_params, transcript_dict, glyph_dict, variant_axes, line_axes, variant_ls, user_tracks, user_track_glyphs, user_lines, user_line_glyphs, title=''):

    project_coords.adjust_coordinates(transcript_dict['exons'], intron_size=plot_params['intron_size'])
    plot_height = plot_params['plot_height']    
    plot = figure(title=title, plot_width=1500, tools='tap,box_zoom,xpan,reset', plot_height=plot_height,min_border=0,#, toolbar_location=None,
               x_range=Range1d(0, transcript_dict['exons'][-1]['compact_end']), y_range=Range1d(0,plot_height), background_fill_color='white')
    
    plot.grid.grid_line_color = None
    plot.toolbar.active_drag = None
    plot.yaxis[0].visible = False

    transcript_ID = transcript_dict['ID'].split(':')[-1]
    
    ### VARIANTS ###

    if variant_ls:
        project_coords.map_point(variant_ls, transcript_dict['exons'])
        try:
            ray_glyph,circle_glyph = add_variant_glyph(plot_params, variant_params, transcript_ID, plot, variant_ls)
        except:
            print(add_variant_glyph(plot_params, variant_params, transcript_ID, plot, variant_ls))

        if ray_glyph and circle_glyph:
            tooltips_variants = [('Position (compact)', '@x'), ('Position (chr)', '@pos'), ('Severity', '@sev')]

            if variant_params['variant_format'].lower().strip('.') == 'vcf': 
                tooltips_variants.extend([('Allele count', '@allele_counts'), ('Allele number', '@allele_numbers'), ('Allele frequency', '@allele_frequencies'), ('Ref', '@ref'), ('Alt', '@alt')])
                tooltips_variants.extend([(info_field, '@'+info_field) for info_field in variant_params['info_annotations']]) #add user defined INFO fields to hover box annotation
                if variant_params['vep']: tooltips_variants.extend([(vep_field, '@'+vep_field) for vep_field in variant_params['vep']['vep_fields']]) #add user defined VEP field to hover box annotation

            plot.add_tools(HoverTool(tooltips=tooltips_variants, renderers=[circle_glyph,ray_glyph], point_policy='follow_mouse', attachment='below'))

        glyph_dict['Variant'].extend([ray_glyph,circle_glyph])

        if variant_params['add_variant_axis']:
            variant_ls = [v for v in variant_ls if v['compact_pos'] >= 0]

            def log10(f): return np.log10(f) if f > 0 else 0

            allele_counts = [v['allele_count'] for v in variant_ls]
            allele_numbers = [v['allele_number'] for v in variant_ls]
            allele_frequencies = [v['allele_frequency'] for v in variant_ls]

            variant_axes['count_linear'].append(add_variant_axis(plot_params, variant_params, plot, 'Allele count', allele_counts, visible=(variant_params['default_y_axis'] == 'AC' and variant_params['default_y_axis_scale'] == 'linear')))
            variant_axes['count_log'].append(add_variant_axis(plot_params, variant_params, plot, 'Log(Allele count)', [log10(c) for c in allele_counts], visible=(variant_params['default_y_axis'] == 'AC' and variant_params['default_y_axis_scale'] == 'log')))
            variant_axes['frequency_linear'].append(add_variant_axis(plot_params, variant_params, plot, 'Allele frequency', allele_frequencies, visible=(variant_params['default_y_axis'] == 'AF' and variant_params['default_y_axis_scale'] == 'linear')))
            variant_axes['frequency_log'].append(add_variant_axis(plot_params, variant_params, plot, 'Log(Allele frequency)', [log10(f) for f in allele_frequencies], visible=(variant_params['default_y_axis'] == 'AF' and variant_params['default_y_axis_scale'] == 'log'))) #don't include 0s since they will just be plotted as min

    else: variant_params['add_variant_axis'] = False

    tooltips_features = [('Type','@feat_type'), ('Start (compact)', '@adj_start'), ('End (compact)', '@adj_end'), 
                      ('Start (chr)', '@true_start'), ('End (chr)', '@true_end'), ('Length', '@true_len')]

    ### INTRONS ###
    introns = project_coords.get_introns_from_exons(transcript_dict['exons'])
    introns_glyph = add_intron_glyph(plot_params, plot, introns)
    plot.add_tools(HoverTool(tooltips=tooltips_features, renderers=[introns_glyph], point_policy='follow_mouse', attachment='below'))    
    
    ### EXONS ###
    exons_glyph,arrow_glyph = add_exon_glyph(plot_params, plot, transcript_dict['exons'], transcript_dict['direction'])
    plot.add_tools(HoverTool(tooltips=tooltips_features, renderers=[exons_glyph], point_policy='follow_mouse', attachment='below'))
    
    glyph_dict['Direction'].append(arrow_glyph)
    glyph_dict['exon'].append(exons_glyph)

    ### UTRs ###
    project_coords.map_box(transcript_dict['UTRs'], transcript_dict['exons'])
    UTR_glyph = add_UTR_glyph(plot_params, plot, transcript_dict['UTRs'])
    plot.add_tools(HoverTool(tooltips=tooltips_features, renderers=[UTR_glyph], point_policy='follow_mouse', attachment='below'))    
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
        plot.add_tools(HoverTool(tooltips=tooltips_tracks, renderers=[track_glyph], point_policy='follow_mouse', attachment='below'))
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

    for sev in variant_params['variant_severity_colors']: variant_params['variant_severity_colors'][sev] = get_color(variant_params['variant_severity_colors'][sev])

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
    try: #if info_annotations empty or nonexistent, set to empty list
        if not variant_params['info_annotations']: variant_params['info_annotations'] = []
    except: variant_params['info_annotations'] = []

    try: #if vep empty or nonexistent, set to empty list
        if not variant_params['vep']: variant_params['vep'] = []
    except: variant_params['vep'] = []

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
    variant_ls = get_variants(variant_params, transcripts, gene_feature.start, gene_feature.end) if variant_params['plot_variants'] else []
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
    
    for idx,ID in enumerate(transcript_IDs):
        title = 'gene={}; transcript={}/{}'.format(plot_params['gene_name'], ID, transcripts[ID]['ID'])
        if transcripts[ID]['direction']: title += ' ({})'.format(transcripts[ID]['direction'])
        plot,glyph_dict = plot_transcript(plot_params, variant_params, user_track_params, user_line_params, transcripts[ID], glyph_dict, variant_axes, line_axes, variant_ls, user_tracks, user_track_glyphs, user_lines, user_line_glyphs, title=title)
        plot_ls.append(plot)

    legend = [add_legend(user_line_params)] if user_line_params else []

    if output_format == 'html': #only add widgets for HTML  

        empty_plot = figure(plot_height=1500,outline_line_color=None) #for HTML, add white space at bottom so hover boxes are not cut off 
        empty_plot.line(x=[0], y=[0]) #avoid empty plot warning
        empty_plot.yaxis.visible = empty_plot.xaxis.visible = empty_plot.grid.visible = False
        plot_ls.append(empty_plot) 

        from bokeh.plotting import output_file, save

        add_exon_zoom(plot_ls,glyph_dict)
        checkbox = add_checkbox(plot_ls,line_axes,glyph_dict,plot_params, variant_params)
        user_tracks_checkbox = add_user_tracks_checkbox(plot_ls,{**variant_axes,**line_axes},user_track_glyphs,glyph_dict['Direction'],plot_params)
        
        user_line_checkboxes=[] #one checkbox per user line, so that they can be lined up with sliders
        for axis in user_line_params:
            user_line_checkbox = add_user_lines_checkbox(plot_ls, line_axes[axis], user_line_glyphs[axis], axis)
            user_line_checkboxes.append(user_line_checkbox)
        
        if variant_params['plot_variants']:
            div_type,radio_group_type,div_scale,radio_group_scale = add_linear_log_scale(variant_params, variant_axes, glyph_dict)
        
        sliders = []
        for axis_name in user_line_params:
            if user_line_params[axis_name]['smoothing_slider']:
                fill_area_ls = [user_line_params[axis_name]['lines'][line]['fill_area'] for line in user_line_params[axis_name]['lines']]*len(plot_ls)
                sliders.append(add_smoothing_slider(user_line_glyphs[axis_name], fill_area_ls, title='{} smoothing'.format(axis_name)))
            else:
                sliders.append(None)
        
        if variant_params['plot_variants'] and variant_params['add_variant_axis']:
            grid1 = [[checkbox, user_tracks_checkbox],[div_type, div_scale],[radio_group_type, radio_group_scale]] 
        else: grid1 = [[checkbox, user_tracks_checkbox]]
        
        lines = list(zip(user_line_checkboxes,sliders))
        for tup in lines: grid1.append(tup)
        grid1.append(legend)
        grid = gridplot(grid1, toolbar_location=None)
        output_file(output)
        save(column([grid]+plot_ls))

    elif output_format == 'png':
        from bokeh.io import export_png
        for p in plot_ls: p.toolbar_location = None
        export_png(column(legend + plot_ls), filename=output)

    elif output_format == 'svg':
        from bokeh.io import export_svg
        for p in plot_ls: p.output_backend = 'svg'
        export_svg(column(legend + plot_ls), filename=output)

main = transcriptionary
