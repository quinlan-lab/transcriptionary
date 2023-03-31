import gffutils
from . import project_coords

def gff_to_db(gff_path,output_name):
    try:
        db = gffutils.FeatureDB(gff_path, keep_order=True)
    except:
        gffutils.create_db(gff_path, dbfn=output_name, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        #db = gffutils.FeatureDB(gff_path, keep_order=True)
        db = gffutils.FeatureDB(output_name, keep_order=True)
    return db

def get_gene_feature(gff_db, gene_name):
    #get ID from name
    for f in gff_db.all_features(featuretype='gene'):
        try:
            if gene_name in f['Name'][0]: return f
        except: continue
    raise ValueError('No gene with name {}'.format(gene_name))

def get_transcript_dict(plot_params, gff_db, gene_feature, transcript_IDs):
    chr_num = gene_feature.seqid
    transcript_dict = {}

    #get all exons and flatten
    if 'flattened-exons' in transcript_IDs or transcript_IDs=='all':
        exons_all = gff_db.children(gene_feature, featuretype='exon')
        exons_all = [{'start': e.start, 'end': e.end} for e in exons_all]
        exons_all = project_coords.flatten_exons(exons_all)
        transcript_dict['flattened-exons'] = dict(ID='flattened-exons', chr_num=chr_num, exons=exons_all, direction='', UTRs=[])    
        
    def get_UTRs(transcript):
        if not plot_params['plot_UTRs']: return []
        UTRs = []
        for f in ['five_prime_UTR', 'three_prime_UTR']:
            try:
                UTR = list(gff_db.children(transcript, featuretype=f))
                if len(UTR) > 1: 
                    print('Warning: transcript {} has {} {} annotations; using first annotation.'.format(tname, len(UTR), featuretype))
                UTR_dict = {'featuretype': f, 'start': UTR[0].start, 'end': UTR[0].end}
                UTRs.append(UTR_dict)
            except: continue 
        return UTRs
                
    #indiv mRNAs
    transcripts = list(gff_db.children(gene_feature, featuretype='mRNA')) + list(gff_db.children(gene_feature, featuretype='transcript'))
    possible_transcripts = [t['Name'][0] for t in transcripts] + [t['ID'][0] for t in transcripts] + [t['ID'][0].split(':')[-1] for t in transcripts] #user can access by Name (ex. KCNQ2-201), ID (ex. transcript:ENST00000344425), or the ID after the colon (ex. ENST00000344425)
    if transcript_IDs != 'all':
        for ID in transcript_IDs:
            if ID not in possible_transcripts and ID != 'flattened-exons': print('No such transcript {}; skipping'.format(ID))
        transcripts = [t for t in transcripts if t['Name'][0]  in transcript_IDs or t['ID'][0] in transcript_IDs or t['ID'][0].split(':')[-1] in transcript_IDs]

    for t in transcripts:
        exons = list(gff_db.children(t, featuretype='exon'))
        exon_coords = [{'start': e.start, 'end': e.end, 'compact_start': -1, 'compact_end': -1} for e in exons]
        direction = t.strand if plot_params['plot_direction'] else ''
        tname = t['Name'][0]
        UTRs = get_UTRs(t)
        transcript_dict[tname] = dict(ID=t['ID'][0], chr_num=chr_num, exons=exon_coords, direction=direction, UTRs=UTRs)
    return transcript_dict
