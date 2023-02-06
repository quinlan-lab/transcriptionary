import csv
import tabix

import itertools as it
try:
    izip = it.izip
except AttributeError:
    izip = zip
    basestring = str

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

        from cyvcf2 import VCF
        from geneimpacts import VEP, Effect

        vep_keys = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info'.split('|')

        vcf = VCF(vcf_path)

        variant_ls = []

        for v in vcf('{}:{}-{}'.format(seqid,start,end)):

            vep = v.INFO['vep']
            impact_strings = vep.split(",")

            impacts = [VEP(impact_string, keys=vep_keys) for impact_string in impact_strings]
            top = Effect.top_severity(impacts)

            try: ann = top[0]
            except: ann = top

            variant_ls.append(dict(pos=v.POS, compact_pos=-1, ref=v.REF, alt=v.ALT,
                    annotation=ann.gene, severity=ann.impact_severity,
                    allele_count=v.INFO.get('AC'),
                    allele_number=v.INFO.get('AN'),
                    allele_frequency=v.INFO.get('AF')))

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
    boxes = []
    db = user_track_params[track_name]['db']
    seqid = user_track_params[track_name]['seqid']
    for s in list(db.region(seqid=seqid, featuretype='exon')):
        box_dict = dict(ID=s['gene_id'][0], start=s.start, end=s.end, compact_start=-1,
                        compact_end=-1, strand=s.strand)  # more than one gene name?
        boxes.append(box_dict)

    return boxes
