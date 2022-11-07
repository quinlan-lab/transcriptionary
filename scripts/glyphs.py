import numpy as np
from bokeh.models import ColumnDataSource, Rect, Segment, Circle, MultiPolygons
from colors import lighten_hex_color

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
        xs_dict.append([{'exterior': adjust_arrow_coords(i, [-20, -20, 0, 0], direction), 'holes': []}]) #rectangle
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
                                      fill_alpha='fill_alpha', line_alpha=0)
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
    hover_color = lighten_hex_color(color, 50)

    glyph = Rect(x='x', y='y', width='w', height='h', height_units='screen', fill_color=color, fill_alpha=fill_alpha,
                 line_color=color)
    hover_glyph = Rect(x='x', y='y', width='w', height='h', height_units='screen', fill_color=hover_color,
                       fill_alpha=fill_alpha, line_color=hover_color)

    return plot.add_glyph(source, glyph, hover_glyph=hover_glyph)


def add_UTR_glyph(plot_params, plot, UTRs, fill_alpha=0.4):
    features_original = [(UTR['start'], UTR['end']) for UTR in UTRs]
    features_compact = [(UTR['compact_start'], UTR['compact_end']) for UTR in UTRs]
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

def add_variant_glyph(plot_params, variant_params, plot, variant_ls, line_width=2):
    variant_ls = [v for v in variant_ls if v['compact_pos'] >= 0]
    N = len(variant_ls)
    if N == 0: 
        variant_params['add_variant_axis'] = False
        return None, None, [0], [0]
    allele_counts = [v['allele_count'] for v in variant_ls]
    allele_numbers = [v['allele_number'] for v in variant_ls]
    allele_frequencies = [v['allele_frequency'] for v in variant_ls]

    if sum(allele_frequencies) < 0: #TODO case where freq data but not cts
        variant_params['add_variant_axis'] = False
    else: variant_params['add_variant_axis'] = True

    x = [v['compact_pos'] for v in variant_ls]
    y0s = [plot_params['y0'] - plot_params['exon_height'] / 2] * N

    def get_y1(fn, ls):
        ls_fn = [fn(x) for x in ls]
        for idx, x in enumerate(ls_fn):
            if x < 0: raise ValueError('Lollipop height cannot be negative')
        y1_circle = [c * (plot_params['plot_height'] - plot_params['y0'] - variant_params[
            'min_lollipop_height'] - variant_params['lollipop_radius'] - line_width - 2) / max([c for c in ls_fn]) + plot_params['y0'] +
                     variant_params['min_lollipop_height'] for c in ls_fn]
        y1_segment = [y - variant_params['lollipop_radius'] for y in y1_circle]
        return y1_circle, y1_segment

    r = [variant_params['lollipop_radius']] * N
    colors = [v['color'] for v in variant_ls]

    hover_colors = [lighten_hex_color(c, 40) for c in colors]

    if variant_params['add_variant_axis']:
        y1_circle,y1_segment = get_y1(lambda x: x, allele_counts)
    else:
        y1_circle,y1_segment = get_y1(lambda x: x/2, [1 for v in variant_ls])
        y1_circle = [y1/2 for y1 in y1_circle] #if lollipop heights don't hold meaning, set them to 1/2 max height
        y1_segment = [y1/2 for y1 in y1_segment]

    cds_dict = dict(x=x, r=r, y0=y0s, y1_circle=y1_circle, y1_segment=y1_segment, 
                                   allele_counts=allele_counts,
                                   allele_numbers=allele_numbers,
                                   allele_frequencies=allele_frequencies,
                                   pos=[v['pos'] for v in variant_ls],
                                   ref=[v['ref'] for v in variant_ls],
                                   alt=[v['alt'] for v in variant_ls],
                                   ann=[v['annotation'] for v in variant_ls],
                                   sev=[v['severity'] for v in variant_ls],
                                   colors=colors, hover_colors=hover_colors, line_alpha=[1 for v in variant_ls])

    if variant_params['add_variant_axis']:
        y1_ci_li_ct, y1_sg_li_ct = get_y1(lambda x: x, allele_counts)
        y1_ci_lg_ct, y1_sg_lg_ct = get_y1(lambda x: np.log10(x) if x > 0 else 0, allele_counts)
        y1_ci_li_fr, y1_sg_li_fr = get_y1(lambda x: x, allele_frequencies)
        y1_ci_lg_fr, y1_sg_lg_fr = get_y1(lambda x: -np.log10(x) if x > 0 else 0, allele_frequencies)

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
                                      fill_alpha=0, line_color='hover_colors', line_width=line_width + 0.5,
                                      line_alpha='line_alpha')
    circle_glyph.nonselection_glyph = None
    return segment_glyph, circle_glyph, allele_counts, allele_frequencies


def add_track_glyph(plot, tracks, height, y_pos):
    tracks = [di for di in tracks if di['compact_start']>=0]
    tracks = sorted(tracks, key=lambda di: di['end']-di['start'], reverse=True) #make sure smaller tracks are in front of larger tracks
    tracks_compact = [(di['compact_start'], di['compact_end']) for di in tracks]
    tracks_original = [(di['start'], di['end']) for di in tracks]
    track_names = [di['ID'] for di in tracks]
    colors = [di['color'] for di in tracks]

    x = [center_feature(f) for f in tracks_compact]
    y = [y_pos for f in tracks_compact]
    w = [end - start for (start, end) in tracks_compact]
    h = [height for f in tracks_compact]
    hover_colors = [lighten_hex_color(color, 40) for color in colors]

    ### METADATA ###
    adj_start = [start for (start, end) in tracks_compact]
    adj_end = [end for (start, end) in tracks_compact]
    true_start = [start for (start, end) in tracks_original]
    true_end = [end for (start, end) in tracks_original]
    true_len = [end - start + 1 for (start, end) in tracks_original]  # +1 for 1 indexed coords

    source = ColumnDataSource(dict(x=x, y=y, w=w, h=h, colors=colors, hover_colors=hover_colors,
                                   track_names=track_names,
                                   adj_start=adj_start, adj_end=adj_end,
                                   true_len=true_len, true_start=true_start, true_end=true_end))

    track_glyph = plot.rect(source=source, x='x', y='y', width='w', height='h', height_units='screen',
                             fill_color='colors', line_color='colors')
    hover_glyph = Rect(x='x', y='y', width='hover_w', height='h', height_units='screen', fill_color='hover_colors',
                       line_color='hover_colors')
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