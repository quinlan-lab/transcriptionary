import numpy as np
from bokeh.models import ColumnDataSource, Rect, Segment, Circle
from .colors import  color_variants,lighten_hex_color
import math

def flatten(ls): return sum(ls, [])

def center_feature(feature):
    center = feature[0] + (feature[1] - feature[0])/2
    return center

def add_exon_glyph(plot_params, plot, exons, direction):
    features_original = [(exon['start'], exon['end']) for exon in exons]
    features_compact = [(exon['compact_start'], exon['compact_end']) for exon in exons]

    x = [center_feature(f) for f in features_compact]
    y = [plot_params['transcript_height'] for f in features_compact]
    w = [end - start for (start, end) in features_compact]
    h = [plot_params['exon_height'] - 4 for f in features_compact]

    ### METADATA ###
    feat_type = ['Exon' for f in features_compact]
    adj_start = [start for (start, end) in features_compact]
    adj_end = [end for (start, end) in features_compact]
    true_start = [start for (start, end) in features_original]
    true_end = [end for (start, end) in features_original]
    true_len = [end - start + 1 for (start, end) in features_original]  # +1 for 1 indexed coords

    exon_source = ColumnDataSource(dict(x=x, y=y, w=w, h=h,
                                        feat_type=feat_type,
                                        adj_start=adj_start, adj_end=adj_end,
                                        true_len=true_len, true_start=true_start, true_end=true_end))

    color = plot_params['glyph_colors']['exon']
    hover_color = lighten_hex_color(color, 20)

    exon_glyph = Rect(x='x', y='y', width='w', height='h', height_units='screen', fill_color=color, line_color=color)
    hover_glyph = Rect(x='x', y='y', width='hover_w', height='hover_h', height_units='screen', fill_color=hover_color,
                       line_color=hover_color)

    if direction == '':
        exon_glyph = plot.add_glyph(exon_source, exon_glyph, hover_glyph=hover_glyph)
        return exon_glyph, None

    ### ADD ARROWS ###
    def adjust_arrow_coords(i, coords, direction):
        ret = [j * (adj_end[i] - adj_start[i]) / 100 for j in coords]

        transcript_length = max([exon['compact_end'] for exon in exons]) + plot_params['intron_size']
        max_range = 0.01 * transcript_length
        range = max(ret) - min(ret) 
        if range > max_range:
            ret = [r * max_range/range for r in ret]

        if direction == '-': ret = [-r for r in ret]
        return [r + x[i] for r in ret]

    xs_dict = []
    for i in range(len(features_compact)):
        xs_dict.append([{'exterior': adjust_arrow_coords(i, [-20, -20, 1, 1], direction), 'holes': []}]) #rectangle
        xs_dict.append([{'exterior': adjust_arrow_coords(i, [0, 0, 20], direction), 'holes': []}]) #triangle

    ys_dict = []
    for i in range(len(features_compact)):
        ys_dict.append([{'exterior': [plot_params['transcript_height'] - plot_params['exon_height'] / 5,
                                      plot_params['transcript_height'] + plot_params['exon_height'] / 5,
                                      plot_params['transcript_height'] + plot_params['exon_height'] / 5,
                                      plot_params['transcript_height'] - plot_params['exon_height'] / 5],
                         'holes': []}])  # rectangle
        ys_dict.append([{'exterior': [plot_params['transcript_height'] - plot_params['exon_height'] / 2,
                                      plot_params['transcript_height'] + plot_params['exon_height'] / 2,
                                      plot_params['transcript_height']], 
                        'holes': []}])  # triangle

    xs = [[[p['exterior'], *p['holes']] for p in mp] for mp in xs_dict]

    ys = [[[p['exterior'], *p['holes']] for p in mp] for mp in ys_dict]

    color = plot_params['glyph_colors']['arrow']

    arrow_source = ColumnDataSource(dict(x=xs, y=ys, fill_alpha=[1 for i in features_compact] * 2,
                                         direction=[direction for i in features_compact] * 2))
    arrow_glyph = plot.multi_polygons(source=arrow_source, xs='x', ys='y', fill_color=color, line_color=color,
                                      fill_alpha='fill_alpha')
    arrow_glyph.level = 'overlay'
    arrow_glyph.nonselection_glyph = None

    exon_glyph = plot.add_glyph(exon_source, exon_glyph, hover_glyph=hover_glyph, nonselection_glyph=None)

    return exon_glyph, arrow_glyph

def add_intron_glyph(plot_params, plot, introns, fill_alpha=1, width=14):
    features_original = [(intron['start'], intron['end']) for intron in introns]
    features_compact = [(intron['compact_start'], intron['compact_end']) for intron in introns]
    x = [center_feature(f) for f in features_compact]
    y = [plot_params['transcript_height'] for f in features_compact]  # prev 10
    w = [plot_params['intron_size'] * 2 for f in features_compact]  # avoid white spaces
    h = [plot_params['exon_height'] / 2 for f in features_compact]

    ### METADATA ###
    feat_type = ['Intron' for f in features_compact]
    adj_start = [start for (start, end) in features_compact]
    adj_end = [end for (start, end) in features_compact]
    true_start = [start for (start, end) in features_original]
    true_end = [end for (start, end) in features_original]
    true_len = [end - start + 1 for (start, end) in features_original]  # +1 for 1 indexed coords

    source = ColumnDataSource(dict(x=x, y=y, w=w, h=h,
                                   feat_type=feat_type,
                                   adj_start=adj_start, adj_end=adj_end,
                                   true_len=true_len, true_start=true_start, true_end=true_end))

    color = plot_params['glyph_colors']['intron']
    hover_color = lighten_hex_color(color, 30)

    glyph = Rect(x='x', y='y', width='w', height='h', height_units='screen', fill_color=color, fill_alpha=fill_alpha,
                 line_color=color)
    hover_glyph = Rect(x='x', y='y', width='w', height='h', height_units='screen', fill_color=hover_color,
                       fill_alpha=fill_alpha, line_color=hover_color)

    return plot.add_glyph(source, glyph, hover_glyph=hover_glyph, nonselection_glyph=None)


def add_UTR_glyph(plot_params, plot, UTRs, fill_alpha=0.4):
    
    features_original = flatten([[(UTR['start'], UTR['end'])]*len(UTR['compact_start']) for UTR in UTRs]) #duplicate original coordinates for each exon spanned by box
    features_compact = flatten([zip(UTR['compact_start'], UTR['compact_end']) for UTR in UTRs])
    feat_type = [UTR['featuretype'] for UTR in UTRs]
    x = [center_feature(f) for f in features_compact]
    y = [plot_params['transcript_height'] for f in features_compact]
    w = [end - start for (start, end) in features_compact]
    h = [plot_params['exon_height'] - 4 for f in features_compact]

    ### METADATA ###
    adj_start = [start for (start, end) in features_compact]
    adj_end = [end for (start, end) in features_compact]
    true_start = [start for (start, end) in features_original]
    true_end = [end for (start, end) in features_original]
    true_len = [end - start + 1 for (start, end) in features_original]  # +1 for 1 indexed coords

    source = ColumnDataSource(dict(x=x, y=y, w=w, h=h,
                                   feat_type=feat_type,
                                   adj_start=adj_start, adj_end=adj_end,
                                   true_len=true_len, true_start=true_start, true_end=true_end))

    color = plot_params['glyph_colors']['UTR']
    hover_color = lighten_hex_color(color, 40)

    glyph = Rect(x='x', y='y', width='w', height='h', height_units='screen', fill_color=color, fill_alpha=fill_alpha,
                 line_color=color)
    hover_glyph = Rect(x='x', y='y', width='w', height='h', height_units='screen', fill_color=hover_color,
                       fill_alpha=fill_alpha, line_color=hover_color)

    return plot.add_glyph(source, glyph, hover_glyph=hover_glyph)

def add_variant_glyph(plot_params, variant_params, variant_set, transcript_ID, plot, line_width=2): #called for each transcript for each variant set
    variant_ls = [v for v in variant_params[variant_set]['variant_ls'] if v['compact_pos'] >= 0]
    N = len(variant_ls)
    if N == 0: 
        return None, None
    
    x = [v['compact_pos'] for v in variant_ls]
    y0s = [plot_params['y0'] - plot_params['exon_height'] / 2] * N

    def get_y1(ys, all_ys, yaxis, yaxis_scale):

        if yaxis == 'AC' and yaxis_scale == 'linear':
            def round_up_to_half(n): #round float up to 0.5
                return math.ceil(n * 2)/2

            y_max_order = math.floor(math.log10(max(all_ys)))

            min_y_axis = 0
            max_y_axis = round_up_to_half(max(all_ys)*10**(-y_max_order)) * 10**y_max_order

        elif yaxis == 'AF' and yaxis_scale == 'linear':
            min_y_axis = 0
            max_y_axis = 1.0

        elif yaxis_scale == 'log':
            ys = [np.log10(y) if y > 0 else 0 for y in ys]
            all_ys = [np.log10(y) if y > 0 else 0 for y in all_ys]
            min_y_axis = math.floor(min(all_ys))
            max_y_axis = math.ceil(max(all_ys))

            if min_y_axis == max_y_axis:
                if min_y_axis == max_y_axis == 0:
                    min_y_axis = 0
                    max_y_axis = 1
                else:
                    min_y_axis = 0
                    max_y_axis = 2*max_y_axis

            if yaxis == 'AF': # if log AF, convert negative values to positive coordinates
                ys = [y - min_y_axis if y != 0 else 0 for y in ys]

        else: # if lollipop heights are not meaningful
            min_y_axis = 0
            max_y_axis = 1

        y1_circle = [y * (plot_params['plot_height'] - plot_params['y0'] - plot_params['min_lollipop_height'] - plot_params['lollipop_radius'] - line_width - 2) / (max_y_axis - min_y_axis) + plot_params['y0'] + plot_params['min_lollipop_height'] for y in ys]
        y1_segment = [y - plot_params['lollipop_radius'] for y in y1_circle]

        return y1_circle, y1_segment

    r = [plot_params['lollipop_radius']] * N
    color_variants(plot_params, variant_params, variant_ls, transcript_ID)
    colors = [v['color'] for v in variant_ls]
    hover_colors = [lighten_hex_color(c, 40) for c in colors]

    allele_counts = [v['allele_count'] for v in variant_ls]
    allele_frequencies = [v['allele_frequency'] for v in variant_ls]

    all_vars = flatten([variant_params[var_set]['variant_ls'] for var_set in variant_params]) #will need list of all vars to get_y1s since we want every variant set to be plotted on the same scale
    all_vars = [v for v in all_vars if v['compact_pos'] >= 0]
    all_allele_counts = [v['allele_count'] for v in all_vars]
    all_allele_frequencies = [v['allele_frequency'] for v in all_vars]

    #set original lollipop ys based on default_y_axis (allele count or allele freq) and default_y_axis_scale (log or linear)
    if plot_params['default_y_axis'] == 'AC': 
        ys = allele_counts
        all_ys = all_allele_counts
    elif plot_params['default_y_axis'] == 'AF': 
        ys = allele_frequencies
        all_ys = all_allele_frequencies
    
    if plot_params['add_variant_axis']:
        y1_circle,y1_segment = get_y1(ys, all_ys, plot_params['default_y_axis'], plot_params['default_y_axis_scale'])

    else:
        #if lollipop heights don't hold meaning, set them to 1/2 max height
        y1_circle,y1_segment = get_y1([1 for v in variant_ls], all_ys, '', '')
        y1_circle = [y1/2 for y1 in y1_circle] 
        y1_segment = [y1/2 - line_width for y1 in y1_segment]

    cds_dict = dict(x=x, r=r, y0=y0s, y1_circle=y1_circle, y1_segment=y1_segment,
                                variant_set=[v['variant_set'] for v in variant_ls],
                                pos=[v['pos'] for v in variant_ls], sev=[v[transcript_ID + '_severity'] for v in variant_ls],
                                colors=colors, hover_colors=hover_colors, line_alpha=[1 for _ in variant_ls])
    
    try: cds_dict['ref']=[v['ref'] for v in variant_ls]
    except: pass
    try: cds_dict['alt']=[v['alt'] for v in variant_ls]
    except: pass

    if plot_params['add_variant_axis']:
        allele_counts = [v['allele_count'] for v in variant_ls]
        allele_frequencies = [v['allele_frequency'] for v in variant_ls]

        cds_dict['allele_counts']=allele_counts
        cds_dict['allele_frequencies']=allele_frequencies

    for info_field in variant_params[variant_set]['info_annotations']: cds_dict[info_field] = [v[info_field] for v in variant_ls]

    for vep_field in variant_params[variant_set]['vep']['vep_fields']: cds_dict[vep_field] = [v[transcript_ID + '_' + vep_field] for v in variant_ls]

    if plot_params['add_variant_axis']:

        y1_ci_li_ct, y1_sg_li_ct = get_y1(allele_counts, all_allele_counts, 'AC', 'linear')
        y1_ci_lg_ct, y1_sg_lg_ct = get_y1(allele_counts, all_allele_counts, 'AC', 'log')
        y1_ci_li_fr, y1_sg_li_fr = get_y1(allele_frequencies, all_allele_frequencies, 'AF', 'linear')
        y1_ci_lg_fr, y1_sg_lg_fr = get_y1(allele_frequencies, all_allele_frequencies, 'AF', 'log')

        variant_axis_dict = dict(
            y1_ci_li_ct=y1_ci_li_ct, y1_sg_li_ct=y1_sg_li_ct,
            y1_ci_lg_ct=y1_ci_lg_ct, y1_sg_lg_ct=y1_sg_lg_ct,
            y1_ci_li_fr=y1_ci_li_fr, y1_sg_li_fr=y1_sg_li_fr,
            y1_ci_lg_fr=y1_ci_lg_fr, y1_sg_lg_fr=y1_sg_lg_fr)
        cds_dict.update(variant_axis_dict)

    source = ColumnDataSource(cds_dict)

    segment_glyph = plot.segment(source=source, x0='x', y0='y0', x1='x', y1='y1_segment', line_color='colors',
                                 line_width=line_width, line_alpha='line_alpha')
    segment_glyph.hover_glyph = Segment(x0='x', y0='y0', x1='x', y1='y1_segment', line_color='hover_colors',
                                        line_width=line_width, line_alpha='line_alpha')
    segment_glyph.nonselection_glyph = None
    circle_glyph = plot.circle('x', 'y1_circle', source=source, radius='r', radius_units='screen', fill_color='white',
                               fill_alpha=0, line_color='colors', line_width=line_width + 0.5, line_alpha='line_alpha')
    circle_glyph.hover_glyph = Circle(x='x', y='y1_circle', radius='r', radius_units='screen', fill_color='white',
                                      fill_alpha=0, line_color='hover_colors', line_width=line_width + 0.5, line_alpha='line_alpha')
    circle_glyph.nonselection_glyph = None
    return segment_glyph, circle_glyph

def add_track_glyph(user_track_params, track_name, plot, tracks, height, y_pos):
    
    tracks = [di for di in tracks if len(di['compact_start'])>0]
    if not tracks: 
        tracks = [{'ID': 'empty track', 'start': 0, 'end': 0, 'compact_start': [0], 'compact_end': [0], 'color': '#111111', 'strand': ''}]
        for field in user_track_params[track_name]['annotate_with']: tracks[0][field] = ''
    tracks = sorted(tracks, key=lambda di: di['end']-di['start'], reverse=True) #make sure smaller tracks are in front of larger tracks
    tracks_compact = flatten([list(zip(di['compact_start'], di['compact_end'])) for di in tracks])
    color_by = flatten([[di['ID']]*len(di['compact_start']) for di in tracks]) #duplicate color_by for each exon spanned by box
    colors = flatten([[di['color']]*len(di['compact_start']) for di in tracks]) #duplicate color for each exon spanned by box

    x = [center_feature(f) for f in tracks_compact]
    y = [y_pos for f in tracks_compact]
    w = [end - start for (start, end) in tracks_compact]
    h = [height for f in tracks_compact]
    hover_colors = [lighten_hex_color(color, 40) for color in colors]

    ### METADATA ###
    cds_dict = dict(x=x, y=y, w=w, h=h, colors=colors, hover_colors=hover_colors, color_by=color_by,
        adj_start=[coord[0] for coord in tracks_compact], 
        adj_end=[coord[1] for coord in tracks_compact],
        true_start=flatten([[di['start']]*len(di['compact_start']) for di in tracks]), #duplicate original coordinates for each exon spanned by box
        true_end=flatten([[di['end']]*len(di['compact_start']) for di in tracks]), 
        true_len=flatten([[di['end'] - di['start'] + 1]*len(di['compact_start']) for di in tracks]))

    for field in user_track_params[track_name]['annotate_with']: cds_dict[field] = [di[field] for di in tracks]

    source = ColumnDataSource(cds_dict)

    track_glyph = plot.rect(source=source, x='x', y='y', width='w', height='h', height_units='screen', fill_color='colors', line_color='colors')
    hover_glyph = Rect(x='x', y='y', width='hover_w', height='h', height_units='screen', fill_color='hover_colors', line_color='hover_colors')
    track_glyph.hover_glyph = hover_glyph
    track_glyph.nonselection_glyph = None
    return track_glyph

def add_multi_line_glyph(plot_params, plot, xs_ls, ys_ls_unsmoothed, max_y=-1, k=1, y0=0, p=False, fill_area=False,
                         y_axis_name=None, line_color='black', line_width=2, line_alpha=1):
    def smooth_ls(ls, k=1):
        ls_copy = ls
        for i in range(0, len(ls_copy), k):
            window_avg = np.mean(ls_copy[i:i + k])
            for j in range(i, i + k):
                if j < len(ls_copy): ls_copy[j] = float(window_avg)
        return ls_copy

    if max_y < 0: max_y = max([y for ys in ys_ls_unsmoothed for y in ys])

    ys_ls_unsmoothed = [[y * (plot_params['plot_height'] - y0) / max_y + y0 for y in ys] for ys in ys_ls_unsmoothed]

    ys_ls_smoothed = [smooth_ls(ys, k=k) for ys in ys_ls_unsmoothed]

    line_alpha = [line_alpha for xs in xs_ls]

    if fill_area:
        xs_ls = [[xs[0]] + xs + [xs[-1]] for xs in xs_ls]
        ys_ls_smoothed = [[y0] + ys + [y0] for ys in ys_ls_smoothed]
        ys_ls_unsmoothed = [[y0] + ys + [y0] for ys in ys_ls_unsmoothed]
        source = ColumnDataSource(dict(x=xs_ls, y=ys_ls_smoothed, y_unsmoothed=ys_ls_unsmoothed, line_alpha=line_alpha))
        multi_line_glyph = plot.patches(source=source, xs='x', ys='y', fill_color=line_color, fill_alpha='line_alpha',
                                        line_alpha=0)
    else:
        source = ColumnDataSource(dict(x=xs_ls, y=ys_ls_smoothed, y_unsmoothed=ys_ls_unsmoothed, line_alpha=line_alpha))
        multi_line_glyph = plot.multi_line(source=source, xs='x', ys='y', line_color=line_color, line_width=line_width,
                                           line_alpha='line_alpha')

    multi_line_glyph.nonselection_glyph = None
    return multi_line_glyph
