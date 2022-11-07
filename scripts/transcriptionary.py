from process_gene_gff import gff_to_db, get_gene_feature, get_transcript_dict
from get_coords import get_variants,get_line,get_track
from map_coords import map_line
from colors import color_boxes, color_variants
from axes import add_user_axis,add_variant_axis
from glyphs import add_intron_glyph, add_exon_glyph, add_variant_glyph, add_UTR_glyph, add_track_glyph, add_multi_line_glyph
from widget_callbacks import add_checkbox,add_user_tracks_checkbox,add_user_lines_checkbox,add_smoothing_slider,add_legend,add_linear_log_scale,add_exon_zoom
import project_coords
import numpy as np
import argparse
from bokeh.plotting import figure, output_file, save
from bokeh.layouts import column, row, gridplot
from bokeh.models import ColumnDataSource, Range1d, HoverTool, LabelSet
import yaml
from yaml.loader import SafeLoader

def constraint_view_plot(plot_params, variant_params, user_line_params, transcript_dict, glyph_dict, axes, variant_ls, user_tracks, user_track_glyphs, user_lines, user_line_glyphs, title=''):

    project_coords.adjust_coordinates(transcript_dict['exons'], intron_size=plot_params['intron_size'])
    plot_height = plot_params['plot_height']    
    plot = figure(title=title, plot_width=1500, tools='tap,box_zoom,xpan,reset', plot_height=plot_height,min_border=0,#, toolbar_location=None,
               x_range=Range1d(0, transcript_dict['exons'][-1]['compact_end']), y_range=Range1d(0,plot_height), background_fill_color='white')
    
    plot.grid.grid_line_color = None
    plot.toolbar.active_drag = None
    plot.yaxis[0].visible = False
    
    ### VARIANTS ###
    if variant_ls:
        project_coords.map_point(variant_ls, transcript_dict['exons'])
        ray_glyph,circle_glyph,allele_counts,allele_frequencies = add_variant_glyph(plot_params, variant_params, plot, variant_ls)
        if ray_glyph and circle_glyph:
            tooltips_variations = [('Position (compact)', '@x'), ('Position (chr)', '@pos'), 
                                   ('Allele count', '@allele_counts'), ('Allele number', '@allele_numbers'), ('Allele frequency', '@allele_frequencies'), 
                                   ('Change', '@ref > @alt'), ('VEP Annotation', '@ann'), ('Severity', '@sev')]
            plot.add_tools(HoverTool(tooltips=tooltips_variations, renderers=[circle_glyph,ray_glyph], point_policy='follow_mouse', attachment='below'))
        glyph_dict['Variant'].extend([ray_glyph,circle_glyph])
        
        if variant_params['add_variant_axis']:
            axes['count'].append(add_variant_axis(plot_params, variant_params, plot, 'Allele count', allele_counts, visible=True))
            axes['allele_frequency'].append(add_variant_axis(plot_params, variant_params, plot, 'Allele frequency',allele_frequencies, visible=False))

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

        tooltips_tracks = [('Name','@track_names'), ('Start (adjusted)', '@adj_start'), ('End (adjusted)', '@adj_end'), 
                          ('Start (true)', '@true_start'), ('End (true)', '@true_end'), ('Length', '@true_len')]
        
        y = ((h*1.5)*(len(user_tracks) - idx - 1)+(h))
        track_glyph = add_track_glyph(plot, user_tracks[track_name], h*0.9, y)
        plot.add_tools(HoverTool(tooltips=tooltips_tracks, renderers=[track_glyph], point_policy='follow_mouse', attachment='below'))
        user_track_glyphs[track_name].append(track_glyph)
        cs = ColumnDataSource(dict(x=[0], y=[(h*1.5)*(len(user_tracks) - idx - 1)], text=[list(user_tracks.keys())[idx]]))
        label = LabelSet(source=cs, x='x', y='y', text='text',text_font_size='{}px'.format(plot_params['track_height']),
                    text_align='left')
        plot.add_layout(label)
        user_track_glyphs[track_name].append(label)
        
    ### USER LINES ###
    for idx,axis_name in enumerate(user_line_params):
        all_xs = []
        all_ys = []
        for line in user_line_params[axis_name]['lines']:
                xs_ls,ys_ls = map_line(user_lines[axis_name][line], transcript_dict['exons'])
                all_xs.append(xs_ls)
                all_ys.append(ys_ls)
        y_max = max([item for sublist in list(np.array(all_ys,dtype=object).flat) for item in sublist])
        
        for idx,line in enumerate(user_line_params[axis_name]['lines']):
            line_params=user_line_params[axis_name]['lines'][line]
            line_glyph = add_multi_line_glyph(plot_params, plot, all_xs[idx], all_ys[idx], max_y=y_max, y0=plot_params['y0'], fill_area=line_params['fill_area'], line_color=line_params['color'], line_alpha=line_params['alpha'])
            line_glyph.level = 'underlay'
            user_line_glyphs[axis_name].append(line_glyph)
            
        add_user_axis(plot, plot_params, user_line_params, axis_name, y_max, plot_params['y0'], plot_params['plot_height'], y_min=0, num_ticks=3, axis_position='right', visible=True)

    return plot,glyph_dict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help='must be .yaml; sample config file can be found at test/test.yaml')
    parser.add_argument('transcripts', help='transcripts to plot; can be [... str], "all", or "transcript_names" to print list of transcripts')
    parser.add_argument('-o', '--output', required=False, help='output filepath')
    args = parser.parse_args()

    ### CONFIG ###
    config_file = args.config_file
    # TODO what to do if config file doesn't have required params?
    if config_file[-5:] != '.yaml': raise ValueError('configuration file must be .yaml')
    with open(config_file) as f:
        params = list(yaml.load_all(f, Loader=SafeLoader))
        plot_params = params[0]
        variant_params = params[1]
        user_track_params = params[2]
        user_line_params = params[3]
        plot_params['transcript_height'] = 20 + len(user_track_params) * plot_params['track_height'] * 1.5
        plot_params['y0'] = plot_params['transcript_height'] + plot_params['exon_height'] / 2  # "0" for line plots

    ### COLORS ###
    with open('default_colors/named_colors.yaml') as f:
        named_colors = list(yaml.load_all(f, Loader=SafeLoader))[0]

    ### GLYPH COLORS ###
    for glyph_type in plot_params['glyph_colors']:
        color = plot_params['glyph_colors'][glyph_type]
        if color[0] != "#": plot_params['glyph_colors'][glyph_type] = named_colors[color]

    ### TRACKS ###
    for track_name in user_track_params:
        track_db = gff_to_db(user_track_params[track_name]['gtf_path'], user_track_params[track_name]['gtf_path'] + '.db')
        user_track_params[track_name]['db'] = track_db
    with open('default_colors/palettes.yaml') as f:
        palettes = list(yaml.load_all(f, Loader=SafeLoader))[0]
    track_colors = [named_colors[c] for c in palettes[plot_params['track_palette']]]
    plot_params['track_colors'] = track_colors

    ### TRANSCRIPTS ###
    transcript_IDs = args.transcripts
    #if transcript_IDs != 'all':
    if transcript_IDs not in ['all', 'transcript_names']:
        transcript_IDs = transcript_IDs.strip('[').strip(']').split(',')
        transcript_IDs = [t.strip('\'') for t in transcript_IDs]

    ### OUTPUT ###
    if args.output: output = args.output + '.html' if args.output[-5:] != '.html' else args.output
    else: output = 'plot.html'
    return plot_params,variant_params,user_track_params,user_line_params,transcript_IDs,output

def transcriptionary():
    plot_params, variant_params, user_track_params, user_line_params, transcript_IDs, output = parse_args()
    # 'test.db' to make sure .db does not accidentally get overwritten during development    
    gff_db = gff_to_db(plot_params['gff_path'],plot_params['gff_path']+'test.db')
    
    gene_feature = get_gene_feature(gff_db, plot_params['gene_name'])

    if transcript_IDs == 'transcript_names':
        print('db1')
        transcripts = get_transcript_dict(plot_params, gff_db, gene_feature, 'all')
        print(list(transcripts.keys()))
        exit()

    transcripts = get_transcript_dict(plot_params, gff_db, gene_feature, transcript_IDs)
    transcript_IDs = list(transcripts.keys()) #if 'all', transcript_IDs will become list of transcript names; if nonexistent IDs they are removed
    
    variant_ls = get_variants(plot_params['variant_path'], gene_feature.start, gene_feature.end, variant_params['seqid']) if variant_params['plot_variants'] else []
    color_variants(plot_params, variant_params, variant_ls)
    user_tracks = {track_name: get_track(user_track_params, track_name) for track_name in user_track_params}
    for track_name in user_tracks: color_boxes(plot_params, user_track_params, track_name, user_tracks[track_name])
    user_lines = {axis_name:{} for axis_name in user_line_params}
    for axis_name in user_line_params:
        for line_name in user_line_params[axis_name]['lines']:
            user_lines[axis_name][line_name] = get_line(user_line_params[axis_name]['lines'][line_name]['filepath'])
    
    plot_ls = []
    glyph_dict = dict(exon=[],UTRs=[],Variant=[],Direction=[])
    user_track_glyphs = {track_name:[] for track_name in user_track_params}
    user_line_glyphs = {line_name:[] for line_name in user_line_params}
    axes = dict(count=[],allele_frequency=[])
    for line_name in user_line_params:
        axes[line_name] = []
    
    for idx,ID in enumerate(transcript_IDs):
        title = 'gene={}; transcript={}/{}'.format(plot_params['gene_name'], ID, transcripts[ID]['ID'])
        if transcripts[ID]['direction']: title += ' ({})'.format(transcripts[ID]['direction'])
        plot,glyph_dict = constraint_view_plot(plot_params, variant_params, user_line_params, transcripts[ID], glyph_dict, axes, variant_ls, user_tracks, user_track_glyphs, user_lines, user_line_glyphs, title=title)
        plot_ls.append(plot)
    
    legend = add_legend(user_line_params)
    add_exon_zoom(plot_ls,glyph_dict)
    checkbox = add_checkbox(plot_ls,axes,glyph_dict,plot_params, variant_params)
    user_tracks_checkbox = add_user_tracks_checkbox(plot_ls,axes,user_track_glyphs,glyph_dict['Direction'],plot_params)
    
    user_line_checkboxes=[] #one checkbox per user line, so that they can be lined up with sliders
    for axis in user_line_params:
        user_line_checkbox = add_user_lines_checkbox(plot_ls, axes[axis], user_line_glyphs[axis], axis)
        user_line_checkboxes.append(user_line_checkbox)
    
    if variant_params['plot_variants']:
        div_type,radio_group_type,div_scale,radio_group_scale = add_linear_log_scale(axes, glyph_dict)
    
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
    grid1.append([legend])
    grid = gridplot(grid1, toolbar_location=None)
    
    output_file(output)
    save(column([grid]+plot_ls))

transcriptionary()