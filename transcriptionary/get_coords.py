import csv
import pandas as pd
from .process_gene_gff import gff_to_db

def get_variants(variant_params, transcripts, start, end):
    filepath = variant_params['filepath']
    seqid = variant_params['chrom']

    transcript_IDs = [transcripts[key]['ID'].split(':')[-1] for key in transcripts]

    def get_variants_vcf(vcf_path):

        from cyvcf2 import VCF
        from geneimpacts import VEP, Effect

        vcf = VCF(vcf_path)
        if variant_params['vep']: vep_keys = vcf.get_header_type(variant_params['vep']['field_name'])['Description'].strip('"').split('Format: ')[-1].split('|')

        variant_ls = []

        for v in vcf('{}:{}-{}'.format(seqid,start,end)):

            di_variant = dict(pos=v.POS, compact_pos=-1, ref=v.REF, alt=v.ALT,            
                    allele_count=v.INFO.get('AC'),
                    allele_number=v.INFO.get('AN'),
                    allele_frequency=v.INFO.get('AF'))

            for info_field in variant_params['info_annotations']: di_variant[info_field] = v.INFO.get(info_field)

            if variant_params['vep']:

                vep = v.INFO[variant_params['vep']['field_name']]
                impact_strings = vep.split(",")

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
                    if len(transcript_impact_string) > 1: print('Warning: more than one VEP impact string with transcript ID {}; using first impact string'.format(transcript_ID))
                    if len(transcript_impact_string) == 0: 
                        di_variant[transcript_ID + '_severity'] = 'NONE'
                    else: 
                        transcript_impact_string = transcript_impact_string[0]
                        di_variant[transcript_ID + '_severity'] = VEP(transcript_impact_string, keys=vep_keys).impact_severity

            else:
                for transcript_ID in transcript_IDs: di_variant[transcript_ID + '_severity'] = 'NONE'

            variant_ls.append(di_variant)

        return variant_ls

    def get_variants_bed(bed_path):
        variant_ls = []

        df = pd.read_csv(bed_path, names=variant_params['header'], sep='\t') 
        for idx,row in df.iterrows():
            if row.iloc[0] != variant_params['chrom']: continue
            di_variant = dict(pos=row.iloc[1], compact_start=-1)
            for transcript_ID in transcript_IDs: di_variant[transcript_ID + '_severity'] = 'NONE'
            variant_ls.append(di_variant)

        return variant_ls
    
    form = variant_params['variant_format'].lower().strip('.')
    if form == 'vcf':
        variant_ls = get_variants_vcf(filepath)

    elif form == 'bed':
        variant_ls = get_variants_bed(filepath)
        
    else: raise ValueError('Invalid filepath: {}'.format(filepath))

    return variant_ls

def get_line(user_line_params,axis_name,line_name):

    filepath = user_line_params[axis_name]['lines'][line_name]['filepath']
    form = user_line_params[axis_name]['lines'][line_name]['format'].lower().strip('.')
    seqid = user_line_params[axis_name]['lines'][line_name]['chrom'].lower().strip('.')

    fi = open(filepath,'r')

    if form == 'bedgraph' or form == 'tsv':
        line_coords = [line.strip().split('\t') for line in fi.readlines()]
    elif form == 'csv':
        line_coords = [line.strip().split(',') for line in fi.readlines()]
    else:
        raise ValueError('Invalid file format: ' + filepath)

    if form == 'bedgraph':
        line_coords = [{'chrom': chrom,'start': int(start), 'end': int(end), 'y': float(y)} for (chrom, start, end, y) in line_coords]

    elif form == 'tsv' or form == 'csv':
        line_coords = [{'chrom': chrom, 'pos': float(x), 'y': float(y)} for (chrom, x, y) in line_coords]
    
    line_coords = [coord for coord in line_coords if coord['chrom'] == seqid]

    fi.close()
    return line_coords

def get_track(user_track_params, track_name):
    form = user_track_params[track_name]['format'].lower()

    boxes = []

    if form in ['gtf', 'gff']:
        db = gff_to_db(user_track_params[track_name]['filepath'], user_track_params[track_name]['filepath'] + '.db')
        seqid = user_track_params[track_name]['chrom']
        for s in list(db.region(seqid=seqid, featuretype='exon')):
            box_dict = dict(ID=s[user_track_params[track_name]['color_by']][0], start=s.start, end=s.end, compact_start=-1, compact_end=-1, strand=s.strand)
            for field in user_track_params[track_name]['annotate_with']:
                box_dict[field] = s[field]
            boxes.append(box_dict)

        user_track_params[track_name]['annotate_with'].append('strand') #automatically annotate gtf tracks with strand

    elif form == 'bed': # does not support headers
        df = pd.read_csv(user_track_params[track_name]['filepath'], names=user_track_params[track_name]['header'], sep='\t')
        for idx,row in df.iterrows():
            if row.iloc[0] != user_track_params[track_name]['chrom']: continue
            box_dict = dict(ID=row[user_track_params[track_name]['color_by']], start=row.iloc[1], end=row.iloc[2], compact_start=-1, compact_end=-1)
            for field in user_track_params[track_name]['annotate_with']:
                box_dict[field] = row[field]
            boxes.append(box_dict)

    return boxes
