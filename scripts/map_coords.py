from project_coords import map_box,map_point

#x coords must be ints
def map_line(line_coords,exons,chr_num=-1):
    xs_ls = []
    ys_ls = []
    
    if 'start' in line_coords[0].keys():
        map_box(line_coords, exons)
        for co in line_coords:
            if co['compact_start'] < 0: continue
            xs = list(range(co['compact_start'],co['compact_end']+1))
            ys = []
            for x in xs:
                ys.append(co['y'])    
            xs_ls.append(xs)
            ys_ls.append(ys)
    else:   
        map_point(line_coords, exons)

        for exon in exons:
            exon_pts = [co for co in line_coords if exon['compact_start'] <= co['compact_pos'] <= exon['compact_end']]
            xs_ls.append([co['compact_pos'] for co in exon_pts])
            ys_ls.append([co['y'] for co in exon_pts])    
    return xs_ls,ys_ls