import pandas as pd
from .process_gene_gff import gff_to_db
from cyvcf2 import VCF
from geneimpacts import VEP, Effect

def get_variants(plot_params, variant_params, variant_set, transcripts, start, end):
    filepath = variant_params[variant_set]['filepath']
    seqid = variant_params[variant_set]['chrom']

    transcript_IDs = [transcripts[key]['ID'].split(':')[-1] for key in transcripts]

    def get_variants_vcf(vcf_path):

        vcf = VCF(vcf_path)

        variant_ls = []

        try:
            vcf.get_header_type('AC')
            vcf.get_header_type('AF')
            variant_params[variant_set]['has_yaxis_info'] = True
        except:
            print('Variant set \'{}\' does not have all both \'AC\' and \'AF\' fields, so lollipop heights will be set to 0'.format(variant_set))
            variant_params[variant_set]['has_yaxis_info'] = False

        def get_allele_field(v, allele_field): #if vcf does not have AC/AN/AF field return 0
            return v.INFO.get(allele_field) if v.INFO.get(allele_field) else 0

        for v in vcf('{}:{}-{}'.format(seqid,start,end)):

            di_variant = dict(pos=v.POS, compact_pos=-1, ref=v.REF, alt=';'.join(v.ALT),            
                    allele_count=get_allele_field(v,'AC'),
                    allele_number=get_allele_field(v,'AN'),
                    allele_frequency=get_allele_field(v,'AF'),
                    variant_set=variant_set)
            
            for info_field in variant_params[variant_set]['info_annotations']: di_variant[info_field] = v.INFO.get(info_field)

            if variant_params[variant_set]['vep']['field_name']:
                
                vep_keys = vcf.get_header_type(variant_params[variant_set]['vep']['field_name'])['Description'].strip('"').split('Format: ')[-1].split('|')
                plot_params['add_variant_severity_checkbox'] = True

                vep = v.INFO[variant_params[variant_set]['vep']['field_name']]
                impact_strings = vep.split(",")

                if variant_params[variant_set]['vep']['annotate_severity_by'] == 'max_severity': #if max_severity, get annotation with highest severity
                        impacts = [VEP(impact_string, keys=vep_keys) for impact_string in impact_strings]
                        top = Effect.top_severity(impacts)

                        try: max_sev_ann = top[0]
                        except: max_sev_ann = top

                #add user defined VEP fields
                for vep_field in variant_params[variant_set]['vep']['vep_fields']:
                    for transcript_ID in transcript_IDs:
                        if variant_params[variant_set]['vep']['annotate_severity_by'] == 'max_severity': impact_string = impact_string = [s for s in impact_strings if max_sev_ann.transcript in s] #if max_severity, use annotation with highest severity
                        elif variant_params[variant_set]['vep']['annotate_severity_by'] == 'transcript_severity': impact_string = [s for s in impact_strings if transcript_ID in s] #elif transcript_severity, use annotation for this transcript
                        if len(impact_string) > 1: print('Warning: more than one VEP impact string with transcript ID {}; using first impact string'.format(transcript_ID))
                        if len(impact_string) == 0: 
                            di_variant[transcript_ID + '_' + vep_field] = 'None'
                        else: 
                            impact_string = impact_string[0]
                            di_variant[transcript_ID + '_' + vep_field] = VEP(impact_string, keys=vep_keys).effects.get(vep_field, None)
                
                #add severity annotation from VEP
                for transcript_ID in transcript_IDs:
                    if variant_params[variant_set]['vep']['annotate_severity_by'] == 'max_severity': impact_string = impact_string = [s for s in impact_strings if max_sev_ann.transcript in s] #if max_severity, use annotation with highest severity
                    elif variant_params[variant_set]['vep']['annotate_severity_by'] == 'transcript_severity': impact_string = [s for s in impact_strings if transcript_ID in s] #elif transcript_severity, use annotation for this transcript
                    if len(impact_string) > 1: print('Warning: more than one VEP impact string with transcript ID {}; using first impact string'.format(transcript_ID))
                    if len(impact_string) == 0: 
                        di_variant[transcript_ID + '_severity'] = 'NONE'
                    else: 
                        impact_string = impact_string[0]
                        di_variant[transcript_ID + '_severity'] = VEP(impact_string, keys=vep_keys).impact_severity

            else:
                for transcript_ID in transcript_IDs: di_variant[transcript_ID + '_severity'] = 'NONE'

            variant_ls.append(di_variant)

        return variant_ls

    def get_variants_bed(bed_path):
        variant_ls = []

        with open(bed_path, 'r') as fi:
            if fi.readline().count('\t') != len(variant_params[variant_set]['header']) - 1:
                print('Warning: number of fields in variant set {} header argument not equal to number of fields in file {} (check for trailing delimiters).'.format(variant_set, bed_path))

        df = pd.read_csv(bed_path, names=variant_params[variant_set]['header'], sep='\t', index_col=None) 
        variant_params[variant_set]['has_yaxis_info'] = False

        for _,row in df.iterrows():
            if str(row.iloc[0]) != str(variant_params[variant_set]['chrom']): continue
            di_variant = dict(pos=row.iloc[1], compact_start=-1, variant_set=variant_set, info_annotations=variant_params[variant_set]['info_annotations'], vep=variant_params[variant_set]['vep'], allele_count=0, allele_frequency=0, allele_number=0)
            if variant_params[variant_set]['consequence_idx']: #0-based
                plot_params['add_variant_severity_checkbox'] = True
                for transcript_ID in transcript_IDs: 
                    di_variant[transcript_ID + '_severity'] = VEP(row.iloc[variant_params[variant_set]['consequence_idx']],['Consequence']).impact_severity
            else:
                for transcript_ID in transcript_IDs:
                    di_variant[transcript_ID + '_severity'] = 'NONE'
            variant_ls.append(di_variant)

        return variant_ls
    
    form = variant_params[variant_set]['format'].lower().strip('.')
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
        for _,row in df.iterrows():
            if row.iloc[0] != user_track_params[track_name]['chrom']: continue
            box_dict = dict(ID=row[user_track_params[track_name]['color_by']], start=row.iloc[1], end=row.iloc[2], compact_start=-1, compact_end=-1)
            for field in user_track_params[track_name]['annotate_with']:
                box_dict[field] = row[field]
            boxes.append(box_dict)

    return boxes
