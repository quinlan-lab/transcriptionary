import numpy as np
import random

def lighten_hex_color(color, delta):
    def limit(x):
        if x > 255: return 255
        if x < 0: return 0
        return x

    if '#' in color: color = color[1:]

    r = hex(limit(int(color[:2], 16) + delta))[2:]
    g = hex(limit(int(color[2:4], 16) + delta))[2:]
    b = hex(limit(int(color[4:6], 16) + delta))[2:]

    new_color = '#{}{}{}'.format(r, g, b)
    return new_color

def color_variants(plot_params, variant_params, variant_ls, transcript_ID):
    for v in variant_ls:
        if v[transcript_ID + '_severity'] != 'NONE': v['color'] = variant_params['variant_severity_colors'][v[transcript_ID + '_severity']]
        else: v['color'] = plot_params['glyph_colors']['variant']
        
def color_boxes(plot_params, user_track_params, track_name, boxes):
    IDs = np.unique([b['ID'] for b in boxes])
    color_dict = {ID: random.choice(np.unique(plot_params['track_colors'])) for ID in IDs}

    try:
        for (ID, color) in user_track_params[track_name]['colors'].items():
            color_dict[ID] = color
    except: pass

    for idx in range(len(boxes)):
        boxes[idx]['color'] = color_dict[boxes[idx]['ID']]