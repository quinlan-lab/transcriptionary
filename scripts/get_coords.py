import csv

def variant_type(annotation):
    if annotation in ['synonymous_variant']:
        return 'LOW'
    if annotation in ['missense_variant', 'inframe_deletion', 'inframe_insertion','stop_lost',]:
        return 'MODERATE'
    if annotation in ['splice_acceptor_variant','frameshift_variant','splice_donor_variant','stop_gained',]:
        return 'HIGH'
    if annotation in ['3_prime_UTR_variant', '5_prime_UTR_variant', 'splice_region_variant','intron_variant',]:
        return 'MODIFIER'
    #print('No classification for variant type {}; classified as "other"'.format(annotation))
    return 'MODIFIER'

def get_variants(filepath, start, end, seqid):
    def get_variants_vcf(vcf_path):
        variant_ls = []
        import tabix #putting import here for now because it doesn't work on local machine
        tb = tabix.open(vcf_path)
        #TODO cast seqid in parse_args
        #for idx,r in enumerate(tb.query(str(seqid),start-1400000,end-1400000)):
        for idx,r in enumerate(tb.query(str(seqid),start,end)):
            #variant_ls.append(dict(pos=int(r[1])+1400000, compact_pos=-1, ref=r[3], alt=r[4], #adding 1.4m bc this vcf is on old build
            variant_ls.append(dict(pos=int(r[1]), compact_pos=-1, ref=r[3], alt=r[4],
                              annotation=r[7].split(';')[-1].split('|')[1], severity=r[7].split(';')[-1].split('|')[2],
                              allele_count=int(r[7].split(';')[0].split('=')[1]),
                              allele_number=int(r[7].split(';')[1].split('=')[1]),
                              allele_frequency=float(r[7].split(';')[2].split('=')[1])))
        return variant_ls

    def get_variants_csv(csv_path):
        variant_ls = []
        with open(csv_path, newline='') as csvfile:
            reader = csv.reader(csvfile)
            next(reader) #skip header
            for row in reader:
                variant_dict = dict(pos=int(row[1]), compact_pos=-1, ref=row[3], alt=row[4], 
                               annotation=row[12], severity=variant_type(row[12]), 
                               allele_count=int(row[16]), allele_number=int(row[17]), allele_frequency=float(row[18]))
                variant_ls.append(variant_dict)
        variant_ls = list({v['pos']:v for v in variant_ls}.values()) #remove duplicate dicts
        return variant_ls

    def get_variants_bed(bed_path):
        import tabix #putting import here for now because it doesn't work on local machine
        variant_ls = []
        tb = tabix.open(bed_path)
        #TODO cast seqid in parse_args
        #for idx,r in enumerate(tb.query(str(seqid),start-1400000,end-1400000)):
        for r in tb.query(str(seqid),start,end):
            variant_ls.append(dict(pos=int(r[1]), compact_pos=-1, ref=r[3], alt=r[4], #adding 1.4m bc this bed is on old build
                             annotation=r[8], severity=variant_type(r[8]),
                             allele_count=-1,
                             allele_number=-1,
                             allele_frequency=-1))
        return variant_ls
    
    if filepath.split('.')[-2] == 'vcf':
        variant_ls = get_variants_vcf(filepath)
    
    elif filepath.split('.')[-1] == 'csv':
        variant_ls = get_variants_csv(filepath)

    #TODO what to do about tabix'd vs non tabix'd?
    elif '.bed' in filepath:
        variant_ls = get_variants_bed(filepath)
        
    else: raise ValueError('Invalid filepath: {}'.format(filepath))

    return variant_ls

def get_line(filepath,chr_num=-1):
    fi = open(filepath,'r')
    line_coords = [line.strip().split('\t') for line in fi.readlines()]
    
    #boxes: <chr>\t<start>\t<end>\t<y>
    if filepath.split('.')[-1] == 'bedgraph' or len(line_coords[0]) == 4:
        line_coords = [{'chrom': int(chrom),'start': int(start), 'end': int(end), 'y': float(y)} for (chrom, start, end, y) in line_coords]
        if chr_num >= 0: coords = [coord for coord in line_coords if coord['chrom']==chr_num]
            
    elif filepath.split('.')[-1] == 'txt':
        if len(line_coords[0]) == 3:
            line_coords = [{'chrom': chrom, 'pos': float(x), 'y': float(y)} for (chrom, x, y) in line_coords]
            if chr_num >= 0: coords = [coord for coord in line_coords if coord['chrom']==chr_num]
        elif len(line_coords[0]) == 2:
            line_coords = [{'chrom': -1, 'pos': float(x), 'y': float(y)} for (x, y) in line_coords]
        else:
            raise ValueError('Invalid file format')

    fi.close()
    return line_coords


def get_track(user_track_params, track_name, strand='-'):
    boxes = []
    db = user_track_params[track_name]['db']
    seqid = user_track_params[track_name]['seqid']
    for s in list(db.region(seqid=seqid, featuretype='exon', strand=strand)):
        box_dict = dict(ID=s['gene_id'][0], start=s.start, end=s.end, compact_start=-1,
                        compact_end=-1)  # more than one gene name?
        boxes.append(box_dict)

    return boxes
