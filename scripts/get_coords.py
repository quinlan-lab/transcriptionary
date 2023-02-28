import csv
import pandas as pd
import tabix
from process_gene_gff import gff_to_db

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

def get_variants(variant_params, transcripts, start, end):
    filepath = variant_params['variant_path']
    seqid = variant_params['seqid']

    transcript_IDs = [transcripts[key]['ID'].split(':')[-1] for key in transcripts]

    def get_variants_vcf(vcf_path):

        from cyvcf2 import VCF
        from geneimpacts import VEP, Effect

        vcf = VCF(vcf_path)
        vep_keys = vcf.get_header_type(variant_params['vep']['field_name'])['Description'].strip('"').split('Format: ')[-1].split('|')

        variant_ls = []

        for v in vcf('{}:{}-{}'.format(seqid,start,end)):

            vep = v.INFO[variant_params['vep']['field_name']]
            impact_strings = vep.split(",")

            di_variant = dict(pos=v.POS, compact_pos=-1, ref=v.REF, alt=v.ALT,            
                    allele_count=v.INFO.get('AC'),
                    allele_number=v.INFO.get('AN'),
                    allele_frequency=v.INFO.get('AF'))

            for info_field in variant_params['info_annotations']: di_variant[info_field] = v.INFO.get(info_field)

            for vep_field in variant_params['vep']['vep_fields']:
                for transcript_ID in transcript_IDs:
                    transcript_impact_string = [s for s in impact_strings if transcript_ID in s]
                    if len(transcript_impact_string) > 1: print('Warning: more than one VEP impact string with transcript ID {}; using first impact string'.format(transcript_ID))
                    if len(transcript_impact_string) == 0: 
                        di_variant[transcript_ID + '_' + vep_field] = 'None'
                    else: 
                        transcript_impact_string = transcript_impact_string[0]
                        di_variant[transcript_ID + '_' + vep_field] = VEP(transcript_impact_string, keys=vep_keys).effects.get(vep_field, 'None')

            for transcript_ID in transcript_IDs:
                transcript_impact_string = [s for s in impact_strings if transcript_ID in s]
                if len(transcript_impact_string) > 1: print(transcript_impact_string)
                if len(transcript_impact_string) == 0: 
                    di_variant[transcript_ID + '_severity'] = 'NONE'
                else: 
                    transcript_impact_string = transcript_impact_string[0]
                    di_variant[transcript_ID + '_severity'] = VEP(transcript_impact_string, keys=vep_keys).impact_severity

            variant_ls.append(di_variant)

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
        variant_ls = []
        tb = tabix.open(bed_path)
        #TODO cast seqid in parse_args
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
    ext = filepath.split('.')[-1]

    if ext == 'bedgraph' or ext == 'tsv':
        line_coords = [line.strip().split('\t') for line in fi.readlines()]
    elif ext == 'csv':
        line_coords = [line.strip().split(',') for line in fi.readlines()]
    else:
        raise ValueError('Invalid file format: ' + filepath)

    if ext == 'bedgraph':
        line_coords = [{'chrom': chrom,'start': int(start), 'end': int(end), 'y': float(y)} for (chrom, start, end, y) in line_coords]
        if chr_num >= 0: coords = [coord for coord in line_coords if coord['chrom']==chr_num]

    elif ext == 'tsv' or ext == 'csv':
        line_coords = [{'chrom': chrom, 'pos': float(x), 'y': float(y)} for (chrom, x, y) in line_coords]
        if chr_num >= 0: coords = [coord for coord in line_coords if coord['chrom']==chr_num]

    fi.close()
    return line_coords

def get_track(user_track_params, track_name):
    filepath = user_track_params[track_name]['filepath'].lower().strip('.')
    form = user_track_params[track_name]['format'].lower()

    boxes = []

    if form in ['gtf', 'gtfdb', 'gff', 'gffdb']:
        
        db = gff_to_db(user_track_params[track_name]['filepath'], user_track_params[track_name]['filepath'] + '.db')
        seqid = user_track_params[track_name]['seqid']
        for s in list(db.region(seqid=seqid, featuretype='exon')):
            # box_dict = dict(ID=s['gene_id'][0], start=s.start, end=s.end, compact_start=-1, compact_end=-1, strand=s.strand)
            box_dict = dict(ID=s[user_track_params[track_name]['color_by']][0], start=s.start, end=s.end, compact_start=-1, compact_end=-1, strand=s.strand)
            for field in user_track_params[track_name]['annotate_with']:
                box_dict[field] = s[field]
            boxes.append(box_dict)

        user_track_params[track_name]['annotate_with'].append('strand') #automatically annotate gtf tracks with strand

    elif form == 'bed': # does not support headers
        df = pd.read_csv(user_track_params[track_name]['filepath'], names=user_track_params[track_name]['header'], sep='\t')
        for idx,row in df.iterrows():
            box_dict = dict(ID=row[user_track_params[track_name]['color_by']], start=row['start'], end=row['end'], compact_start=-1, compact_end=-1)
            for field in user_track_params[track_name]['annotate_with']:
                box_dict[field] = row[field]
            boxes.append(box_dict)

    return boxes