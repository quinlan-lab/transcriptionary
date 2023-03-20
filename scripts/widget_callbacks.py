from bokeh.models import ColumnDataSource,CustomJS,CheckboxGroup,Legend,Slider,Div,RadioGroup
from bokeh.plotting import figure
from axes import format_ticks
import numpy as np

def add_linear_log_scale(variant_params, axes, glyph_dict):
    div_type = Div(text="""Lollipop height:""", width=200, height=15)
    div_scale = Div(text="""Lollipop height scale:""", width=200, height=15)
    
    variant_type_labels = ['Allele count', 'Allele frequency']
    variant_scale_labels = ['Linear', 'Log']
    
    def log10(f): return np.log10(f) if f > 0 else 0

    radio_group_type = RadioGroup(labels=variant_type_labels, active=0 if variant_params['default_y_axis'] == 'AC' else 1, width=200)
    radio_group_scale = RadioGroup(labels=variant_scale_labels, active=0 if variant_params['default_y_axis_scale'] == 'linear' else 1, width=200)

    radio_group_type.js_on_click(CustomJS(args=dict(radio_group_scale=radio_group_scale,axes=axes, 
                                                    variant=glyph_dict['Variant']), 
        code="""        
        function setAxisVis(axis_name, vis){
            for (let i = 0; i < axes[axis_name].length; i++){
                axes[axis_name][i].visible = vis
            }
        }
        
        for (let i = 0; i < variant.length; i++){
            if (cb_obj.active == 0) {
                if (radio_group_scale.active == 0){
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_li_ct']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_li_ct']
                    setAxisVis('count_linear', true)
                    setAxisVis('count_log', false)
                    setAxisVis('frequency_linear', false)
                    setAxisVis('frequency_log', false)
                } else {
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_lg_ct']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_lg_ct']
                    setAxisVis('count_linear', false)
                    setAxisVis('count_log', true)
                    setAxisVis('frequency_linear', false)
                    setAxisVis('frequency_log', false)
                }
            }
            if (cb_obj.active == 1) {
                if (radio_group_scale.active == 0){
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_li_fr']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_li_fr']
                    setAxisVis('count_linear', false)
                    setAxisVis('count_log', false)
                    setAxisVis('frequency_linear', true)
                    setAxisVis('frequency_log', false)
                } else {
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_lg_fr']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_lg_fr']
                    setAxisVis('count_linear', false)
                    setAxisVis('count_log', false)
                    setAxisVis('frequency_linear', false)
                    setAxisVis('frequency_log', true)
                }
            }
            variant[i].data_source.change.emit()        
        }
    """))    
        
    radio_group_scale.js_on_click(CustomJS(args=dict(radio_group_type=radio_group_type,axes=axes, 
                                                     variant=glyph_dict['Variant']), code="""        
        function setAxisVis(axis_name, vis){
            for (let i = 0; i < axes[axis_name].length; i++){
                axes[axis_name][i].visible = vis
            }
        }
        
        for (let i = 0; i < variant.length; i++){
            if (cb_obj.active == 0) {
                if (radio_group_type.active == 0){
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_li_ct']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_li_ct']
                    setAxisVis('count_linear', true)
                    setAxisVis('count_log', false)
                    setAxisVis('frequency_linear', false)
                    setAxisVis('frequency_log', false)
                    
                } else {
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_li_fr']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_li_fr']
                    setAxisVis('count_linear', false)
                    setAxisVis('count_log', false)
                    setAxisVis('frequency_linear', true)
                    setAxisVis('frequency_log', false)
                }
            }
            if (cb_obj.active == 1) {
                if (radio_group_type.active == 0){
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_lg_ct']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_lg_ct']
                    setAxisVis('count_linear', false)
                    setAxisVis('count_log', true)
                    setAxisVis('frequency_linear', false)
                    setAxisVis('frequency_log', false)
                    
                } else {
                    variant[i].data_source.data['y1_circle']=variant[i].data_source.data['y1_ci_lg_fr']
                    variant[i].data_source.data['y1_segment']=variant[i].data_source.data['y1_sg_lg_fr']
                    setAxisVis('count_linear', false)
                    setAxisVis('count_log', false)
                    setAxisVis('frequency_linear', false)
                    setAxisVis('frequency_log', true)
                    
                }
            }
            variant[i].data_source.change.emit()        
        }
    """))     
    
    return div_type,radio_group_type,div_scale,radio_group_scale

def add_checkbox(plot_ls, line_axes, glyph_dict, plot_params, variant_params):
    labels = []
    for tup in zip(['UTRs','Direction','Variant'], [plot_params['plot_UTRs'], plot_params['plot_direction'], variant_params['plot_variants']]):
        if tup[1]: labels.append(tup[0])
    active = list(range(len(labels)))
    checkbox = CheckboxGroup(labels=labels,active=active, width=100)
    
    glyph_dict_checkbox = {label:glyph_dict[label] for label in labels}

    callback = CustomJS(args=dict(plot_ls=plot_ls, line_axes=line_axes, labels=labels, glyph_dict_checkbox=glyph_dict_checkbox,
                                 ), code='''
        function setVis(glyphs, vis){
            for (let i = 0; i < glyphs.length; i++){
                if(glyphs[i] != null){
                    glyphs[i].visible = vis
                }
            }
        }
        
        function setAxisVis(axis_name, vis){
            for (let i = 0; i < axes[axis_name].length; i++){
                line_axes[axis_name][i].visible = vis
            }
        }
        
        for (let i = 0; i < labels.length; i++){
            if (cb_obj.active.includes(i)) {
                setVis(glyph_dict_checkbox[labels[i]], true)
                
            }
            else {setVis(glyph_dict_checkbox[labels[i]], false)}
        }
    ''')

    checkbox.js_on_change('active', callback)
    
    return checkbox 

def add_user_tracks_checkbox(plot_ls,axes,user_track_glyphs,direction_glyphs,plot_params):
    labels = list(user_track_glyphs.keys())
    active = list(range(len(labels)))
    user_tracks_checkbox = CheckboxGroup(labels=labels,active=active, width=100)
    axis_names = [ls[0].axis_label for ls in list(axes.values()) if len(ls)>0]
    user_track_names = list(user_track_glyphs.keys())
    
    ori_y_coords = [user_track_glyphs[track_name][0].data_source.data['y'][0] for track_name in user_track_glyphs]
        
    direction_glyphs = [glyph for glyph in direction_glyphs if glyph]
    s0_arrow = [glyph_ls.data_source for glyph_ls in direction_glyphs]
    callback = CustomJS(args=dict(plot_ls=plot_ls, plot_params=plot_params, axis_names=axis_names,axes=axes,
                                  user_track_glyphs=user_track_glyphs, user_track_names=user_track_names,
                                  ori_y_coords=ori_y_coords, s0_arrow=s0_arrow
                                 ), code='''
        const transcript_h = plot_params['transcript_height']
        const exon_h = plot_params['exon_height']
        const track_h = plot_params['track_height']
            
        function setVis(glyphs, vis){
            for (let i = 0; i < glyphs.length; i++){
                if(glyphs[i] != null){
                    glyphs[i].visible = vis
                }
            }
        }
        
        function setAxesStart(start){
            for (let i = 0; i < axis_names.length; i++){
                for (let j = 0; j < plot_ls.length; j++){
                    try {
                        plot_ls[j].y_range.start = start
                        plot_ls[j].extra_y_ranges[axis_names[i]].start = start
                    }
                    catch (err) {}
                }
            }
        }
        
        //in LabelSet source is called source whereas in GlyphRenderer source is called data_source for no apparent reason
        function get_source(glyph){
            try {glyph.data_source.data; return glyph.data_source}
            catch (err) {return glyph.source}
        }
        
        function setY(track_name,y){
            for (let i = 0; i < user_track_glyphs[track_name].length; i++){
                for (let j = 0; j < get_source(user_track_glyphs[track_name][i]).data['y'].length; j++){
                    if (i%2==0) {get_source(user_track_glyphs[track_name][i]).data['y'][j] = y}
                    else {get_source(user_track_glyphs[track_name][i]).data['y'][j] = y-track_h}    
                }
                get_source(user_track_glyphs[track_name][i]).change.emit()
            }
        }
        
        function adjust_arrow(s0_arrow,start){ 
            
            for (let i = 0; i < s0_arrow.length; i++){
                for (let j = 0; j < s0_arrow[i].data['y'].length/2; j++){
                    //triangle
                    s0_arrow[i].data['y'][j*2+1][0][0][0] = transcript_h-exon_h/2 + start/(track_h*1.5) //start 57
                    s0_arrow[i].data['y'][j*2+1][0][0][1] = transcript_h+exon_h/2 - start/(track_h*1.5)
                    //rectangle
                    s0_arrow[i].data['y'][j*2][0][0][0] = transcript_h-exon_h/5 + start/(track_h*1.5*3)
                    s0_arrow[i].data['y'][j*2][0][0][1] = transcript_h+exon_h/5 - start/(track_h*1.5*3)
                    s0_arrow[i].data['y'][j*2][0][0][2] = transcript_h+exon_h/5 - start/(track_h*1.5*3)
                    s0_arrow[i].data['y'][j*2][0][0][3] = transcript_h-exon_h/5 + start/(track_h*1.5*3)
                }
                s0_arrow[i].change.emit()
            }
        }
        
        var temp = 0
        var y_co = 0
        
        for (let i = 0; i < user_track_names.length; i++){
            if (cb_obj.active.includes(i)) {
                setVis(user_track_glyphs[user_track_names[i]],true)
                setY(user_track_names[i], ori_y_coords[temp])
                y_co = ori_y_coords[temp]
                temp++
            }
            else {setVis(user_track_glyphs[user_track_names[i]],false)}
        }
        var start = 0
        if (y_co==0){start = plot_params['transcript_height']-plot_params['exon_height']/2}
        else {start = y_co-ori_y_coords[ori_y_coords.length-1]}
        setAxesStart(start)
        adjust_arrow(s0_arrow,start)
    ''')

    user_tracks_checkbox.js_on_change('active', callback)
    return user_tracks_checkbox

def add_user_lines_checkbox(plot_ls,axis,user_line_glyphs,label):
    user_lines_checkbox = CheckboxGroup(labels=[label],active=[0], width=100, align=('start','center'))
    callback = CustomJS(args=dict(plot_ls=plot_ls, axes=axis,
                                  user_line_glyphs=user_line_glyphs
                                 ), code='''
        function setVis(glyphs, vis){
            for (let i = 0; i < glyphs.length; i++){
                if(glyphs[i] != null) { glyphs[i].visible = vis }
            }
        }
        
        function setAxisVis(vis){
            for (let i = 0; i < axes.length; i++){
                axes[i].visible = vis
            }
        }
                
        if (cb_obj.active.includes(0)) {                
            setVis(user_line_glyphs,true)
            setAxisVis(true)
        }
        else {
            setVis(user_line_glyphs,false)
            setAxisVis(false)
        }
    ''')

    user_lines_checkbox.js_on_change('active', callback)
    return user_lines_checkbox

def add_smoothing_slider(glyph_ls, fill_area_ls, title=''):
    slider = Slider(start=1, end=100, value=25, step=1, title=title, width=250)
    slider = Slider(start=1, end=100, value=1, step=1, title=title, width=250)
    slider.js_on_change("value", CustomJS(args=dict(glyph_ls=glyph_ls, fill_area_ls=fill_area_ls),code="""
        function smooth(y,k,fill_area) {
            var smoothed_ls = [ ... y]
            if (smoothed_ls.length <= 2) {return smoothed_ls}
            if (fill_area) {smoothed_ls = smoothed_ls.slice(1,smoothed_ls.length-1)}
            
            for (let i = 0; i < smoothed_ls.length; i += k){
                const window_avg = mean(smoothed_ls.slice(i, i+k))
                for (let j = i; j < i + k; j++){
                    if (j < smoothed_ls.length-1) { smoothed_ls[j] = window_avg }
                }
            }
            if (fill_area) {return [y[0], ...smoothed_ls, y[y.length-1]]}
            
            else {return smoothed_ls}
        }
        
        function mean(y) { return y.reduce((i, j) => i + j) / y.length }
        
        for (let i = 0; i < glyph_ls.length; i++){
            for (let j = 0; j < glyph_ls[i].data_source.data['y'].length; j++){
                glyph_ls[i].data_source.data['y'][j] = smooth(glyph_ls[i].data_source.data['y_unsmoothed'][j], cb_obj.value, fill_area_ls[i])
            }
            glyph_ls[i].data_source.change.emit()
        }
        
    """))
    return slider

def add_legend(user_line_params, width=270):
    num_lines = sum([len(user_line_params[axis_name]['lines']) for axis_name in user_line_params])
    height = 20 *num_lines + 10
    legend_plot = figure(plot_height=height,plot_width=width,toolbar_location=None)
    legend_plot.yaxis.visible = legend_plot.xaxis.visible = legend_plot.grid.visible = False
    
    legend_items = []
    for axis_name in user_line_params:
        for line in user_line_params[axis_name]['lines']:        
        #TODO add alpha?
            legend_glyph = legend_plot.line(x=[0,0],y=[0,0],line_color=user_line_params[axis_name]['lines'][line]['color'], legend_label=line,)
    
    legend_plot.legend.padding = 0
    legend_plot.legend.spacing = 0
    legend_plot.legend.margin = 0
    legend_plot.legend.location = (0,0)
    legend_plot.legend.border_line_alpha = 0
    legend_plot.outline_line_alpha = 0
    

    return legend_plot

def add_exon_zoom(plot_ls,glyph_dict):
    for (plot,exon_glyph,arrow_glyph) in zip(plot_ls,glyph_dict['exon'], glyph_dict['Direction']):
        s0_exon = exon_glyph.data_source
        s1_exon = ColumnDataSource(data=exon_glyph.data_source.data) #deep copy is not working in CustomJS so I'm doing it here
        try:
            s0_arrow = arrow_glyph.data_source
            s1_arrow = ColumnDataSource(data=arrow_glyph.data_source.data)
        except: s0_arrow = s1_arrow = None

        callback = CustomJS(args=dict(plot=plot, s0_exon=s0_exon, s1_exon=s1_exon, 
                                  s0_arrow=s0_arrow, s1_arrow=s1_arrow,
                       ),code="""
        cb_obj.active = !cb_obj.active
        const N_exon = s1_exon.data['x'].length
        const idx = cb_obj.indices[0]
        const active_exon_start = s1_exon.data['adj_start'][idx]
        const active_exon_end = s1_exon.data['adj_end'][idx]
        const active_exon_len = s1_exon.data['true_len'][idx]
        const plot_len = plot.x_range.end
        
        function adjust_arrow(s0_arrow, s1_arrow){
            const N_arrow = s1_arrow.data['x'].length/2
            var direc_int = 1
            if (s1_arrow.data['direction'][idx*2] == '-') {direc_int = -1}
            const plot_mdpt = active_exon_start + active_exon_len/2
            s0_arrow.data['x'][idx*2][0][0] = [-3, -3, 0, 0].map(function(x) { return x * direc_int * active_exon_len / 200 + plot_mdpt; })
            s0_arrow.data['x'][idx*2+1][0][0] = [0,0,3].map(function(x) { return x * direc_int * active_exon_len / 200 + plot_mdpt; })
        }
        
        function reset_arrow(s0_arrow, s1_arrow){
            const N_arrow = s1_arrow.data['x'].length
            for (let i = 0; i < N_arrow; i++){
                s0_arrow.data['x'][i][0][0] = s1_arrow.data['x'][i][0][0]
            }
        }
        
        function adjust_x_axis() {
            plot.x_range.start = active_exon_start
            plot.x_range.end = active_exon_end
        }
        
        function reset_x_axis() {
            plot.x_range.start = 0
            plot.x_range.end = s1_exon.data['adj_end'][N_exon-1]
        }
                       
        if(cb_obj.active){
            adjust_x_axis()
            if (s1_arrow != null) {adjust_arrow(s0_arrow, s1_arrow)}
        } else{
            reset_x_axis()
            if (s1_arrow != null) {reset_arrow(s0_arrow, s1_arrow)}
        }
        if (s1_arrow != null) {s0_arrow.change.emit()}
    """)

        s0_exon.selected.js_on_change('indices', callback)