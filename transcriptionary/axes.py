from bokeh.models import LinearAxis,Range1d
import math

def add_axis(plot, plot_params, axis_label, y_max, min_tick, max_tick, tick_precision, tick_scientific_notation, y_min=0, num_ticks=3, axis_position='right', visible=True):
    ticks = [(max_tick-min_tick)/(num_ticks+1)*i+min_tick for i in range(num_ticks+2)]
    ticks = list(map(lambda x: int(x) if int(x)==x else x , ticks)) #have to do this for overrides to work on ticks ending in .0
    ticker_dict = {tick:str((y_max-y_min)/(num_ticks+1)*(idx)+y_min) for idx,tick in enumerate(ticks)}
    #tick formatting
    ticker_dict = format_ticks(ticker_dict, tick_precision, tick_scientific_notation)
    plot.extra_y_ranges[axis_label] = Range1d(0,plot_params['plot_height'])
    axis = LinearAxis(y_range_name=axis_label, axis_label=axis_label, ticker=ticks, major_label_overrides=ticker_dict, visible=visible)
    plot.add_layout(axis, axis_position)
    return axis

def add_user_axis(plot, plot_params, user_line_params, axis_name, y_max, min_tick, max_tick, y_min=0, num_ticks=3, axis_position='right', visible=True):
    axis_label = user_line_params[axis_name]['y_axis_label']
    axis = add_axis(plot, plot_params, axis_label, y_max, min_tick, max_tick, user_line_params[axis_name]['tick_precision'], user_line_params[axis_name]['tick_scientific_notation'], y_min=y_min, num_ticks=num_ticks, axis_position=axis_position, visible=visible)
    return axis

def add_variant_axis(plot_params, variant_params, plot, axis_label, allele_vals, visible=True):

    if axis_label == 'Allele count':
        def round_up_to_half(n): #round float up to 0.5
            return math.ceil(n * 2)/2

        y_max_order = math.floor(math.log10(max(allele_vals)))

        y_min = 0
        y_max = round_up_to_half(max(allele_vals)*10**(-y_max_order)) * 10**y_max_order #ex. max allele vals == 1400 --> round 1.4 up to nearest 0.5, then multiply by 1000
        num_ticks = 3

    elif axis_label == 'Allele frequency':
        y_min = 0
        y_max = 1.0
        num_ticks = 3

    elif axis_label in ['Log(Allele count)', 'Log(Allele frequency)']:
        y_min = math.floor(min(allele_vals))
        y_max = math.ceil(max(allele_vals))
        if y_min != y_max:
            num_ticks = y_max - y_min - 1
        elif y_min == y_max == 0:
            y_min = 0
            y_max = 1
            num_ticks = 0
        else:
            y_min = min([0,2*y_max])
            y_max = max([0,2*y_max])
            num_ticks = 1

    min_tick = plot_params['y0']+variant_params['min_lollipop_height']
    max_tick = plot_params['plot_height']-variant_params['lollipop_radius']-variant_params['lollipop_line_width']
    tick_precision = 2
    tick_scientific_notation = False
        
    axis_position = 'left'
    visible = visible

    axis = add_axis(plot, plot_params, axis_label, y_max, min_tick, max_tick, tick_precision, tick_scientific_notation, y_min=y_min, num_ticks=num_ticks, axis_position=axis_position, visible=visible)
    return axis

def format_ticks(ticker_dict, tick_precision, tick_scientific_notation):
    for key in ticker_dict:
        if tick_scientific_notation: format_precision = '.{}e'.format(tick_precision)
        else: format_precision = '.{}f'.format(tick_precision)
        ticker_dict = {key:format(float(value),format_precision) for (key,value) in ticker_dict.items()}
    for key in ticker_dict:
        if '.' in ticker_dict[key]: ticker_dict[key] = ticker_dict[key].rstrip('0').rstrip('.')

    return ticker_dict
