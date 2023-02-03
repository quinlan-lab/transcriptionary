# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 12:40:23 2021

@author: slulla
"""

def adjust_coordinates(exons, intron_size=10):
    
    '''
    UNIT TEST 1 - intron at 0
    >>> ec1 = [{'start': 100, 'end': 200}, {'start': 250, 'end': 400}, {'start': 700, 'end': 900}]
    >>> adjust_coordinates(ec1)
    >>> for exon in ec1: print(exon)
    {'start': 100, 'end': 200, 'compact_start': 10, 'compact_end': 110}
    {'start': 250, 'end': 400, 'compact_start': 120, 'compact_end': 270}
    {'start': 700, 'end': 900, 'compact_start': 280, 'compact_end': 480}
    
    UNIT TEST 2 - exon at 0
    >>> ec2 = [{'start': 0, 'end': 200}, {'start': 250, 'end': 400}, {'start': 700, 'end': 900}]
    >>> adjust_coordinates(ec2)
    >>> for exon in ec2: print(exon)
    {'start': 0, 'end': 200, 'compact_start': 0, 'compact_end': 200}
    {'start': 250, 'end': 400, 'compact_start': 210, 'compact_end': 360}
    {'start': 700, 'end': 900, 'compact_start': 370, 'compact_end': 570}
    
    UNIT TEST 3 - intron at 0, intron_size=20
    >>> ec3 = [{'start': 100, 'end': 200}, {'start': 250, 'end': 400}, {'start': 700, 'end': 900}]
    >>> adjust_coordinates(ec3, intron_size=20)
    >>> for exon in ec3: print(exon)
    {'start': 100, 'end': 200, 'compact_start': 20, 'compact_end': 120}
    {'start': 250, 'end': 400, 'compact_start': 140, 'compact_end': 290}
    {'start': 700, 'end': 900, 'compact_start': 310, 'compact_end': 510}
    
    UNIT TEST 4 - exon at 0, intron_size=20
    >>> ec4 = [{'start': 0, 'end': 200}, {'start': 250, 'end': 400}, {'start': 700, 'end': 900}]
    >>> adjust_coordinates(ec4, intron_size=20)
    >>> for exon in ec4: print(exon)
    {'start': 0, 'end': 200, 'compact_start': 0, 'compact_end': 200}
    {'start': 250, 'end': 400, 'compact_start': 220, 'compact_end': 370}
    {'start': 700, 'end': 900, 'compact_start': 390, 'compact_end': 590}
    
    UNIT TEST 5 - 1 intron length < intron_size
    >>> ec5 = [{'start': 0, 'end': 200}, {'start': 250, 'end': 400}, {'start': 410, 'end': 610}]
    >>> adjust_coordinates(ec5, intron_size=20)
    >>> for exon in ec5: print(exon)
    {'start': 0, 'end': 200, 'compact_start': 0, 'compact_end': 200}
    {'start': 250, 'end': 400, 'compact_start': 220, 'compact_end': 370}
    {'start': 410, 'end': 610, 'compact_start': 380, 'compact_end': 580}
    
     UNIT TEST 6 - backwards coordinates
     >>> ec6 = [{'start': 900, 'end': 700}, {'start': 400, 'end': 250}, {'start': 200, 'end': 100}]
     >>> adjust_coordinates(ec6)
     >>> for exon in ec6: print(exon)
     {'start': 700, 'end': 900, 'compact_start': 280, 'compact_end': 480}
     {'start': 250, 'end': 400, 'compact_start': 120, 'compact_end': 270}
     {'start': 100, 'end': 200, 'compact_start': 10, 'compact_end': 110}
     
     UNIT TEST 7 - empty input
     >>> ec7 = []
     >>> adjust_coordinates(ec7)
     >>> for exon in ec7: print(exon)
     
     UNIT TEST 8 - overlapping exon coordinates
     >>> ec8 = [{'start': 0, 'end': 200}, {'start': 100, 'end': 300}, {'start': 400, 'end': 500}]
     >>> adjust_coordinates(ec8)
     Traceback (most recent call last):
         ...
     ValueError: Exon coordinates overlap
    '''

    def flip_exon_coords(exons):
        for exon in exons:
            tmp = exon['start']
            exon['start'] = exon['end']
            exon['end'] = tmp
    
   
    if len(exons) == 0: return
    
    #if exon coordinates are backwards, flip them
    if exons[0]['start'] > exons[0]['end']: 
        flip_exon_coords(exons)
    
    exons = sorted(exons, key=lambda exon:exon['start'])
    
    end_last_compact = 0 #end of last exon in compacted coords
    end_last_original = 0 #end of last exon in original coords
    
    for idx,exon in enumerate(exons):
        
        if end_last_original > exon['start']: raise ValueError('Exon coordinates overlap')
        
        if exon['start'] == 0:
            exon['compact_start'] = exon['start']
            exon['compact_end'] = exon['end']
            
            end_last_compact = exon['end']
            end_last_original = exon['end']
            continue
        
        intron_len = exon['start'] - end_last_original #length of intron between this exon and previous
        new_start = end_last_compact+intron_size if intron_len > intron_size else end_last_compact+intron_len
        new_end = new_start + exon['end'] - exon['start'] ###1 based or 0 based coordinate system?###
        exon['compact_start'] = new_start
        exon['compact_end'] = new_end
            
        end_last_compact = new_end
        end_last_original = exon['end']

def flatten_exons(exons):
    
    '''
    UNIT TEST 1 - no overlaps
    >>> e1 = [{'start': 150, 'end': 200}, {'start': 0, 'end': 100}, {'start': 300, 'end': 450}, {'start': 500, 'end': 600}]
    >>> flatten_exons(e1)
    [{'start': 0, 'end': 100}, {'start': 150, 'end': 200}, {'start': 300, 'end': 450}, {'start': 500, 'end': 600}]
    
    UNIT TEST 2 - all overlaps
    >>> e2 = [{'start': 50, 'end': 100}, {'start': 75, 'end': 150}, {'start': 0, 'end': 80}, {'start': 120, 'end': 130}]
    >>> flatten_exons(e2)
    [{'start': 0, 'end': 150}]
    UNIT TEST 3 - start with overlaps
    >>> e3 = [{'start': 10, 'end': 50}, {'start': 40, 'end': 70}, {'start': 60, 'end': 80}, {'start': 100, 'end': 150}, {'start': 160, 'end': 200}, {'start': 300, 'end': 350}]
    >>> flatten_exons(e3)
    [{'start': 10, 'end': 80}, {'start': 100, 'end': 150}, {'start': 160, 'end': 200}, {'start': 300, 'end': 350}]
    
    UNIT TEST 4 - end with overlaps
    >>> e4 = [{'start': 50, 'end': 100}, {'start': 150, 'end': 200}, {'start': 250, 'end': 350}, {'start': 375, 'end': 390}, {'start': 325, 'end': 400}, {'start': 270, 'end': 320}, {'start': 400, 'end': 450}, {'start': 380, 'end': 450}, {'start': 450, 'end': 500}]
    >>> flatten_exons(e4)
    [{'start': 50, 'end': 100}, {'start': 150, 'end': 200}, {'start': 250, 'end': 500}]
    
    UNIT TEST 5 - overlaps only in middle
    >>> e5 = [{'start': 200, 'end': 250}, {'start': 275, 'end': 300}, {'start': 350, 'end': 400}, {'start': 310, 'end': 400}, {'start': 360, 'end': 380}, {'start': 400, 'end': 410}, {'start': 450, 'end': 500}]
    >>> flatten_exons(e5)
    [{'start': 200, 'end': 250}, {'start': 275, 'end': 300}, {'start': 310, 'end': 410}, {'start': 450, 'end': 500}]
    '''
        
    exon_coords = [(exon['start'],exon['end']) for exon in exons]
    
    def largest_range(ls):
        min_start = min([start for (start,end) in ls])
        max_end = max([end for (start,end) in ls])
        return (min_start, max_end)

    exon_coords = sorted(exon_coords, key=lambda x:x[0]) #sort by starting coord
    
    exons_flat = [] #list to return
    overlaps = [exon_coords[0]] #list of overlapping exons
    max_end = exon_coords[0][1] #end of overlapping region

    for i in range(1, len(exon_coords)): 
        if exon_coords[i][0] <= max_end: #if start of this exon is before end of overlapping region:
            overlaps.append(exon_coords[i]) #add this exon to overlapping region
            if exon_coords[i][1] > max_end: max_end = exon_coords[i][1]
        else:
            if len(overlaps) == 1: exons_flat.append(overlaps[0])
            if len(overlaps) > 1: exons_flat.append(largest_range(overlaps))
            overlaps = [exon_coords[i]]
            max_end = exon_coords[i][1]

    exons_flat.append(largest_range(overlaps)) #have to do this one more time at the end
    exons_flat = sorted(exons_flat, key=lambda x:x[0]) #sort by starting coord
    
    exons = []
    for idx,e in enumerate(exons_flat):
        exons.append(dict(start=e[0],end=e[1]))    
    return exons



def map_box(box_coords, exons):
    '''
    UNIT TEST 1 - one box inside each exon
    >>> ec1 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc1 = [{'start': 150, 'end': 175}, {'start': 300, 'end': 375}, {'start': 800, 'end': 825}]
    >>> adjust_coordinates(ec1)
    >>> map_box(bc1, ec1)
    >>> for coord in bc1: print(coord)
    {'start': 150, 'end': 175, 'compact_start': 60, 'compact_end': 85}
    {'start': 300, 'end': 375, 'compact_start': 170, 'compact_end': 245}
    {'start': 800, 'end': 825, 'compact_start': 380, 'compact_end': 405}
    
    UNIT TEST 2 - >= 1 box inside each exon 
    >>> ec2 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc2 = [{'start': 150, 'end': 160}, {'start': 170, 'end': 190}, {'start': 270, 'end': 290}, {'start': 300, 'end': 350}, {'start': 360, 'end': 375}, {'start': 800, 'end': 825}, {'start': 830, 'end': 880}]
    >>> adjust_coordinates(ec2)
    >>> map_box(bc2, ec2)
    >>> for coord in bc2: print(coord)
    {'start': 150, 'end': 160, 'compact_start': 60, 'compact_end': 70}
    {'start': 170, 'end': 190, 'compact_start': 80, 'compact_end': 100}
    {'start': 270, 'end': 290, 'compact_start': 140, 'compact_end': 160}
    {'start': 300, 'end': 350, 'compact_start': 170, 'compact_end': 220}
    {'start': 360, 'end': 375, 'compact_start': 230, 'compact_end': 245}
    {'start': 800, 'end': 825, 'compact_start': 380, 'compact_end': 405}
    {'start': 830, 'end': 880, 'compact_start': 410, 'compact_end': 460}
        
    UNIT TEST 3 - 2 exons without a box
    >>> ec3 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc3 = [{'start': 260, 'end': 300}, {'start': 330, 'end': 370}]
    >>> adjust_coordinates(ec3)
    >>> map_box(bc3, ec3)
    >>> for coord in bc3: print(coord)
    {'start': 260, 'end': 300, 'compact_start': 130, 'compact_end': 170}
    {'start': 330, 'end': 370, 'compact_start': 200, 'compact_end': 240}
    
    
    UNIT TEST 4 - boxes partially in exons
    >>> ec4 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc4 = [{'start': 180, 'end': 210}, {'start': 230, 'end': 300}, {'start': 320, 'end': 420}, {'start': 660, 'end': 800}, {'start': 820, 'end': 850}]
    >>> adjust_coordinates(ec4)
    >>> map_box(bc4, ec4)
    >>> for coord in bc4: print(coord)
    {'start': 180, 'end': 210, 'compact_start': 90, 'compact_end': 110}
    {'start': 230, 'end': 300, 'compact_start': 120, 'compact_end': 170}
    {'start': 320, 'end': 420, 'compact_start': 190, 'compact_end': 270}
    {'start': 660, 'end': 800, 'compact_start': 280, 'compact_end': 380}
    {'start': 820, 'end': 850, 'compact_start': 400, 'compact_end': 430}
    
    UNIT TEST 5 - boxes on exon boundaries
    >>> ec5 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc5 = [{'start': 100, 'end': 200}, {'start': 250, 'end': 300}, {'start': 700, 'end': 725}, {'start': 750, 'end': 900}]
    >>> adjust_coordinates(ec5)
    >>> map_box(bc5, ec5)
    >>> for coord in bc5: print(coord)
    {'start': 100, 'end': 200, 'compact_start': 10, 'compact_end': 110}
    {'start': 250, 'end': 300, 'compact_start': 120, 'compact_end': 170}
    {'start': 700, 'end': 725, 'compact_start': 280, 'compact_end': 305}
    {'start': 750, 'end': 900, 'compact_start': 330, 'compact_end': 480}
    
    UNIT TEST 6 - boxes in introns
    >>> ec6 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc6 = [{'start': 0, 'end': 90}, {'start': 210, 'end': 220}, {'start': 230, 'end': 240}, {'start': 1000, 'end': 1100}]
    >>> adjust_coordinates(ec6)
    >>> map_box(bc6, ec6)
    >>> for coord in bc6: print(coord)
    {'start': 0, 'end': 90, 'compact_start': -1, 'compact_end': -1}
    {'start': 210, 'end': 220, 'compact_start': -1, 'compact_end': -1}
    {'start': 230, 'end': 240, 'compact_start': -1, 'compact_end': -1}
    {'start': 1000, 'end': 1100, 'compact_start': -1, 'compact_end': -1}
    
    UNIT TEST 7 - boxes in introns, starting at last base of previous exon and/or ending at last base of next exon
    >>> ec7 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc7 = [{'start': 50, 'end': 100}, {'start': 200, 'end': 210}, {'start': 220, 'end': 250}, {'start': 400, 'end': 500}, {'start': 400, 'end': 700}, {'start': 900, 'end': 950}]
    >>> adjust_coordinates(ec7)
    >>> map_box(bc7, ec7)
    >>> for coord in bc7: print(coord)
    {'start': 50, 'end': 100, 'compact_start': -1, 'compact_end': -1}
    {'start': 200, 'end': 210, 'compact_start': -1, 'compact_end': -1}
    {'start': 220, 'end': 250, 'compact_start': -1, 'compact_end': -1}
    {'start': 400, 'end': 500, 'compact_start': -1, 'compact_end': -1}
    {'start': 400, 'end': 700, 'compact_start': -1, 'compact_end': -1}
    {'start': 900, 'end': 950, 'compact_start': -1, 'compact_end': -1}
    
    UNIT TEST 8 - start and end coordinates reversed
    >>> ec8 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc8 = [{'start': 210, 'end': 180}, {'start': 420, 'end': 320}, {'start': 800, 'end': 660}, {'start': 850, 'end': 820}]
    >>> adjust_coordinates(ec8)
    >>> map_box(bc8, ec8)
    >>> for coord in bc8: print(coord)
    {'start': 180, 'end': 210, 'compact_start': 90, 'compact_end': 110}
    {'start': 320, 'end': 420, 'compact_start': 190, 'compact_end': 270}
    {'start': 660, 'end': 800, 'compact_start': 280, 'compact_end': 380}
    {'start': 820, 'end': 850, 'compact_start': 400, 'compact_end': 430}
    
    UNIT TEST 9 - overlapping boxes
    >>> ec9 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc9 = [{'start': 125, 'end': 175}, {'start': 140, 'end': 190}, {'start': 200, 'end': 300}, {'start': 270, 'end': 370}, {'start': 700, 'end': 850}, {'start': 820, 'end': 870}, {'start': 800, 'end': 900}]
    >>> adjust_coordinates(ec9)
    >>> map_box(bc9, ec9)
    >>> for coord in bc9: print(coord)
    {'start': 125, 'end': 175, 'compact_start': 35, 'compact_end': 85}
    {'start': 140, 'end': 190, 'compact_start': 50, 'compact_end': 100}
    {'start': 200, 'end': 300, 'compact_start': 120, 'compact_end': 170}
    {'start': 270, 'end': 370, 'compact_start': 140, 'compact_end': 240}
    {'start': 700, 'end': 850, 'compact_start': 280, 'compact_end': 430}
    {'start': 820, 'end': 870, 'compact_start': 400, 'compact_end': 450}
    {'start': 800, 'end': 900, 'compact_start': 380, 'compact_end': 480}
    
    UNIT TEST 10 - boxes out of order
    >>> ec10 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc10 = [{'start': 230, 'end': 300}, {'start': 180, 'end': 210}, {'start': 820, 'end': 850}, {'start': 320, 'end': 420}, {'start': 660, 'end': 800}]
    >>> adjust_coordinates(ec10)
    >>> map_box(bc10, ec10)
    >>> for coord in bc10: print(coord)
    {'start': 230, 'end': 300, 'compact_start': 120, 'compact_end': 170}
    {'start': 180, 'end': 210, 'compact_start': 90, 'compact_end': 110}
    {'start': 820, 'end': 850, 'compact_start': 400, 'compact_end': 430}
    {'start': 320, 'end': 420, 'compact_start': 190, 'compact_end': 270}
    {'start': 660, 'end': 800, 'compact_start': 280, 'compact_end': 380}
    
    UNIT TEST 11 - empty box list
    >>> ec11 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc11 = []
    >>> adjust_coordinates(ec11)
    >>> map_box(bc11, ec11)
    >>> for coord in bc11: print(coord)
    
    UNIT TEST 12 - empty exon list
    >>> ec12 = []
    >>> bc12 = [{'start': 150, 'end': 175}, {'start': 300, 'end': 375}, {'start': 800, 'end': 825}]
    >>> adjust_coordinates(ec12)
    >>> map_box(bc12, ec12)
    >>> for coord in bc12: print(coord)
    {'start': 150, 'end': 175, 'compact_start': -1, 'compact_end': -1}
    {'start': 300, 'end': 375, 'compact_start': -1, 'compact_end': -1}
    {'start': 800, 'end': 825, 'compact_start': -1, 'compact_end': -1}
    
    UNIT TEST 13 - >= 1 box inside each exon with intron_size=20
    >>> ec13 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> bc13 = [{'start': 150, 'end': 160}, {'start': 170, 'end': 190}, {'start': 270, 'end': 290}, {'start': 300, 'end': 350}, {'start': 360, 'end': 375}, {'start': 800, 'end': 825}, {'start': 830, 'end': 880}]
    >>> adjust_coordinates(ec13, intron_size=20)
    >>> map_box(bc13, ec13)
    >>> for coord in bc13: print(coord)
    {'start': 150, 'end': 160, 'compact_start': 70, 'compact_end': 80}
    {'start': 170, 'end': 190, 'compact_start': 90, 'compact_end': 110}
    {'start': 270, 'end': 290, 'compact_start': 160, 'compact_end': 180}
    {'start': 300, 'end': 350, 'compact_start': 190, 'compact_end': 240}
    {'start': 360, 'end': 375, 'compact_start': 250, 'compact_end': 265}
    {'start': 800, 'end': 825, 'compact_start': 410, 'compact_end': 435}
    {'start': 830, 'end': 880, 'compact_start': 440, 'compact_end': 490}
    
    
    '''
    
    for b in box_coords:
        if b['start'] > b['end']:
            tmp = b['start']
            b['start'] = b['end']
            b['end'] = tmp
        b['compact_start'] = -1
        b['compact_end'] = -1
        for exon in exons:
            if exon['end'] <= b['start'] or exon['start'] >= b['end']: continue
            b_start_compact = exon['compact_start'] if b['start'] < exon['start'] else exon['compact_start'] + b['start'] - exon['start']
            b_end_compact = exon['compact_end'] if b['end'] > exon['end'] else exon['compact_end'] - (exon['end'] - b['end'])
            b['compact_start'] = b_start_compact
            b['compact_end'] = b_end_compact
            break

def map_point(point_coords, exons):
    
    '''
    UNIT TEST 1 - 1 point inside each exon
    >>> ec1 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc1 = [{'pos': 150}, {'pos': 300}, {'pos': 800},]
    >>> adjust_coordinates(ec1)
    >>> map_point(pc1, ec1)
    >>> for coord in pc1: print(coord)
    {'pos': 150, 'compact_pos': 60}
    {'pos': 300, 'compact_pos': 170}
    {'pos': 800, 'compact_pos': 380}
    
    UNIT TEST 2 - >=1 point inside each exon
    >>> ec2 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc2 = [{'pos': 150}, {'pos': 175}, {'pos': 300}, {'pos': 750}, {'pos': 800}, {'pos': 850}]
    >>> adjust_coordinates(ec2)
    >>> map_point(pc2, ec2)
    >>> for coord in pc2: print(coord)
    {'pos': 150, 'compact_pos': 60}
    {'pos': 175, 'compact_pos': 85}
    {'pos': 300, 'compact_pos': 170}
    {'pos': 750, 'compact_pos': 330}
    {'pos': 800, 'compact_pos': 380}
    {'pos': 850, 'compact_pos': 430}
    
    UNIT TEST 3 - 2 exons without a point
    >>> ec2 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc2 = [{'pos': 300}, {'pos': 350}]
    >>> adjust_coordinates(ec2)
    >>> map_point(pc2, ec2)
    >>> for coord in pc2: print(coord)
    {'pos': 300, 'compact_pos': 170}
    {'pos': 350, 'compact_pos': 220}
    
    UNIT TEST 4 - points on exon boundaries
    >>> ec4 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc4 = [{'pos': 100}, {'pos': 400}, {'pos': 900},]
    >>> adjust_coordinates(ec4)
    >>> map_point(pc4, ec4)
    >>> for coord in pc4: print(coord)
    {'pos': 100, 'compact_pos': 10}
    {'pos': 400, 'compact_pos': 270}
    {'pos': 900, 'compact_pos': 480}
    
    UNIT TEST 5 - points in introns
    >>> ec5 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc5 = [{'pos': 50}, {'pos': 210}, {'pos': 240}, {'pos': 500}, {'pos': 950},]
    >>> adjust_coordinates(ec5)
    >>> map_point(pc5, ec5)
    >>> for coord in pc5: print(coord)
    {'pos': 50, 'compact_pos': -1}
    {'pos': 210, 'compact_pos': -1}
    {'pos': 240, 'compact_pos': -1}
    {'pos': 500, 'compact_pos': -1}
    {'pos': 950, 'compact_pos': -1}
   
    UNIT TEST 6 - overlapping points
    >>> ec6 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc6 = [{'pos': 150}, {'pos': 150}, {'pos': 240}, {'pos': 240}, {'pos': 750}, {'pos': 750}]
    >>> adjust_coordinates(ec6)
    >>> map_point(pc6, ec6)
    >>> for coord in pc6: print(coord)
    {'pos': 150, 'compact_pos': 60}
    {'pos': 150, 'compact_pos': 60}
    {'pos': 240, 'compact_pos': -1}
    {'pos': 240, 'compact_pos': -1}
    {'pos': 750, 'compact_pos': 330}
    {'pos': 750, 'compact_pos': 330}
   
    UNIT TEST 7 - points out of order
    >>> ec7 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc7 = [{'pos': 800}, {'pos': 175}, {'pos': 300}, {'pos': 150}, {'pos': 850}, {'pos': 750}]
    >>> adjust_coordinates(ec7)
    >>> map_point(pc7, ec7)
    >>> for coord in pc7: print(coord)
    {'pos': 800, 'compact_pos': 380}
    {'pos': 175, 'compact_pos': 85}
    {'pos': 300, 'compact_pos': 170}
    {'pos': 150, 'compact_pos': 60}
    {'pos': 850, 'compact_pos': 430}
    {'pos': 750, 'compact_pos': 330}
    
    UNIT TEST 8 - empty point list
    >>> ec8 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc8 = []
    >>> adjust_coordinates(ec8)
    >>> map_point(pc8, ec8)
    >>> for coord in pc8: print(coord)
    
    UNIT TEST 9 - empty exon list
    >>> ec9 = []
    >>> pc9 = [{'pos': 150}, {'pos': 300}, {'pos': 800},]
    >>> adjust_coordinates(ec9)
    >>> map_point(pc9, ec9)
    >>> for coord in pc9: print(coord)
    {'pos': 150, 'compact_pos': -1}
    {'pos': 300, 'compact_pos': -1}
    {'pos': 800, 'compact_pos': -1}
    
    UNIT TEST 10 - >=1 point inside each exon with intron_size=20
    >>> ec10 = [{'start': 100, 'end': 200},{'start': 250, 'end': 400},{'start': 700, 'end': 900}]
    >>> pc10 = [{'pos': 150}, {'pos': 175}, {'pos': 300}, {'pos': 750}, {'pos': 800}, {'pos': 850}]
    >>> adjust_coordinates(ec10, intron_size=20)
    >>> map_point(pc10, ec10)
    >>> for coord in pc10: print(coord)
    {'pos': 150, 'compact_pos': 70}
    {'pos': 175, 'compact_pos': 95}
    {'pos': 300, 'compact_pos': 190}
    {'pos': 750, 'compact_pos': 360}
    {'pos': 800, 'compact_pos': 410}
    {'pos': 850, 'compact_pos': 460}
    
    '''
    
    for coord in point_coords:
        coord['compact_pos'] = -1
        for exon in exons:
            if coord['pos'] >= exon['start'] and coord['pos'] <= exon['end']:
                coord['compact_pos'] = exon['compact_start'] + coord['pos'] - exon['start']
                break

#x coords must be ints
def map_line(line_coords,exons,seqid):
    line_coords = [di for di in line_coords if di['chrom'] == seqid]
    if len(line_coords) == 0: return [[]],[[]]

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

def get_introns_from_exons(exons):
    introns = []

    flat_exons = [(exon['start'], exon['end']) for exon in exons]
    flat_exons = [item for sublist in flat_exons for item in sublist][:-1]
    if len(flat_exons) % 2 == 1: flat_exons = [0] + flat_exons

    flat_exons_compact = [(exon['compact_start'], exon['compact_end']) for exon in exons]
    flat_exons_compact = [item for sublist in flat_exons_compact for item in sublist][:-1]
    if len(flat_exons_compact) % 2 == 1: flat_exons_compact = [0] + flat_exons_compact
    for i in range(0, len(flat_exons) - 1, 2):
        introns.append({'start': flat_exons[i] + 1, 'end': flat_exons[i + 1] - 1,
                        'compact_start': flat_exons_compact[i] + 1, 'compact_end': flat_exons_compact[i + 1] - 1})
    return introns


if __name__ == "__main__":
    import doctest
    doctest.testmod()