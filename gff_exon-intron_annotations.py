import os
import sys
import numpy as np
import matplotlib.pyplot as plt

gff_file = sys.argv[1]
keywords_file = sys.argv[2]
intron_annot_file = sys.argv[3]
exon_intron_annot_file = sys.argv[4]
intron_cords_file = sys.argv[5]
coordinate_type = sys.argv[6] # relative coordinates or absolute coordinates.
median = int(sys.argv[7])

# Depending on the length you want your sequences to be resized to, you can either choose to enter a length here
# alternatively if you want to use the median length for resizing you can comment this out and use the code below on
# median = round_to_nearest_multiple_of_10(median_size)
# To keep it simple, I have used the variable name median and it depends on user to use the actual median size or
# specify a length they'd like the introns to be resized to.

# CAN BE CONVERTED TO BED FILE EASILY WITH EXCEL AND R - crdnts_to_bed.r ###


#Simple logic here is if the exon spans an intron, it likely is a retained intron event whereas if an intron spans an exon, that's not an actual intron and gets filtered out by intron_filter().
def retained_introns(model):
    model['splicing_event'] = np.nan
    exons = model[model['group'] == 'exon']
    introns = model[model['group'] == 'intron']
    exoncoords = set(map(tuple, exons[['start', 'end']].values))

    for index, row in introns.iterrows():
        for exon_start, exon_end in exoncoords:
            if exon_start < row['start'] and exon_end > row['end']:
                model.loc[index, 'splicing_event'] = 'retained intron'
                exon_to_delete = exons[(exons['start'] == exon_start) & (exons['end'] == exon_end)]
                model.loc[exon_to_delete.index, 'splicing_event'] = 'drop'
                    
    # I am dropping those exons which aren't real exons but two exons flanking a retained intron.
    model = model[model['splicing_event'] != 'drop']
    return model


def exon_cat(gene_model):
    group=gene_model[gene_model['group']=='exon']
    group.sort_values(by=['transcript_id'],inplace=True)
    group['transcript_id'] = group['transcript_id'].astype(str).str.replace('transcript_id=','')
    smallest_exon_start_tr = group.groupby(['transcript_id'])[['start','end']].min().reset_index()
    greatest_exon_end_tr = group.groupby(['transcript_id'])[['start','end']].max().reset_index()
    
    group['splicing_event'] = np.nan
    smallest_exon_rows = group[group[['transcript_id','start', 'end']].apply(tuple, axis=1).isin(smallest_exon_start_tr[['transcript_id','start', 'end']].apply(tuple, axis=1))]
    greatest_exon_rows = group[group[['transcript_id','start', 'end']].apply(tuple, axis=1).isin(greatest_exon_end_tr[['transcript_id','start', 'end']].apply(tuple, axis=1))]

    if '+' in group['strand'].unique():
        group.loc[smallest_exon_rows.index, 'splicing_event'] = 'first_exon'
        group.loc[greatest_exon_rows.index, 'splicing_event'] = 'last_exon'
        
    elif '-' in group['strand'].unique():
        group.loc[smallest_exon_rows.index, 'splicing_event'] = 'last_exon'
        group.loc[greatest_exon_rows.index, 'splicing_event'] = 'first_exon'

    excluded = group[
        (group['splicing_event'] == 'first_exon') |
        (group['splicing_event'] == 'last_exon')
        ][['start','end', 'transcript_id']]
    
    excluded_set = set(map(tuple, excluded.values))

    excluded_tuples = group[['start', 'end', 'transcript_id']].apply(tuple, axis=1).isin(excluded_set)
    filtered_group = group[~(excluded_tuples)]
    constitutive_alltrs_dic={}
    tr_numfilter = filtered_group['transcript_id'].nunique()
    alltrexons = filtered_group.groupby(['start','end'])['transcript_id'].count().reset_index(name='transcript_id_count')
    Const_overall_trs_condition = alltrexons[alltrexons['transcript_id_count'] == tr_numfilter]
    if not Const_overall_trs_condition.empty:
        constitutive_alltrs_tuples = set(Const_overall_trs_condition[['start', 'end']].itertuples(index=False, name=None))
    else:
        constitutive_alltrs_tuples=set()

    #Processing skipped and constitutive separately.
    first_exon_rows = group[group['splicing_event'] == 'first_exon']
    first_exon_tuples = first_exon_rows[['start','end','transcript_id']].apply(tuple, axis=1)
    

    # Grouping transcripts by first exons.
    tr_subgroup_df = first_exon_rows.groupby(['start','end'])['transcript_id'].apply(list).reset_index()
    tr_subgroup = tr_subgroup_df['transcript_id'].to_dict()

    skipped_dic = {}
    constitutive_dic = {}
    # group transcripts into sub-groups based on first exon and then locate the exon under question in the transcripts that have first exons coordinates smaller than exon in quesstion !
    # I am creating transcript subgroups with all exons except the first ones. Hence, I use filtered_group for subgroup.
    for key in tr_subgroup:
        subgroup = group[group['transcript_id'].isin(tr_subgroup[key])]
        tr_num = subgroup['transcript_id'].nunique()
        subgroup_filter =subgroup[~(subgroup['splicing_event']=='first_exon') & ~(subgroup['splicing_event']=='last_exon')]
        ex_counts = subgroup_filter.groupby(['start', 'end'])['transcript_id'].count().reset_index(name='transcript_id_count')
        skipped_condition = ex_counts[ex_counts['transcript_id_count'] < tr_num]
        skipped = list(skipped_condition[['start', 'end']].itertuples(index=False, name=None))
        skipped_dic[key] = skipped
        constitutive_condition = ex_counts[ex_counts['transcript_id_count'] == tr_num]
        constitutive = list((constitutive_condition[['start', 'end']].itertuples(index=False, name=None)))
        constitutive_dic[key]=constitutive
        
    # Next approach will be parsing out the true constitutive and skipped based on first exon coordinates.
    constitutive_true_dic={}
    constitutive_tuples_pre = set(tuple_ for tuples_list in constitutive_dic.values() for tuple_ in tuples_list)
    
    for coords in constitutive_tuples_pre:
        start, end = coords
        if '+' in first_exon_rows['strand'].unique():
            # With this line of code I am checking if the exon in question is an internal exon in any of the transcripts. 
            # Only when I locate transcript block by its first exon and then if the exon in question is internal to it ( start > tup[1] (or first_exon end)) then, I extract the transcript_ids
            subset_condition = group.loc[
                group.apply(
                    lambda row: any(start > tup[1]
                        for tup in first_exon_tuples
                        if ((row['end'] == tup[1]) & (row['splicing_event'] == 'first_exon'))
                    ),axis=1)]
            trs = subset_condition['transcript_id'].unique()
            subset_trs = filtered_group[filtered_group['transcript_id'].isin(trs)]
        
        elif '-' in first_exon_rows['strand'].unique():
            subset_condition = group.loc[
                group.apply(
                    lambda row: any(start < tup[1]
                        for tup in first_exon_tuples
                        if ((row['end'] == tup[1]) & (row['splicing_event'] == 'first_exon'))
                    ),axis=1)]
            trs = subset_condition['transcript_id'].unique()
            subset_trs = filtered_group[filtered_group['transcript_id'].isin(trs)]

        second_ex_counts = subset_trs.groupby(['start', 'end'])['transcript_id'].count().reset_index(name='transcript_id_count')
        tr_num_subset = subset_trs['transcript_id'].nunique()
        
        skipped_sec_condition = second_ex_counts[second_ex_counts['transcript_id_count'] < tr_num_subset][['start','end']]
        if start in list(skipped_sec_condition['start']):
            skipped_dic[start]=[]
            skipped_dic[start].append(coords)

        const_real = second_ex_counts[second_ex_counts['transcript_id_count'] == tr_num_subset][['start','end']]
        if start in list(const_real['start']):
            constitutive_true_dic[start]=[]
            constitutive_true_dic[start].append(coords)
      
    skipped_tuples = set(tuple_ for tuples_list in skipped_dic.values() for tuple_ in tuples_list)
    constitutive_tuples = set(tuple_ for tuples_list in constitutive_true_dic.values() for tuple_ in tuples_list)

    #I am accounting for first and last exons here, they need to be excluded when I am trying to classify exons of other transcripts in which the first exon of trs1 might be a second or third exon.
    tr_num = filtered_group['transcript_id'].nunique()
    ex_1start_more_ends = filtered_group.groupby(['start']).filter(lambda g: g['end'].nunique() > 1)
    ex_1start_more_ends_count = ex_1start_more_ends.groupby(['start'])['transcript_id'].transform('count')
    composite_condition1 = ex_1start_more_ends[(ex_1start_more_ends_count < tr_num)]

    ex_1end_more_starts = filtered_group.groupby(['end']).filter(lambda g: g['start'].nunique() > 1)
    ex_1end_more_starts_count = ex_1end_more_starts.groupby(['end'])['transcript_id'].transform('count')
    composite_condition2 = ex_1end_more_starts[(ex_1end_more_starts_count < tr_num)]

    if '+' in group['strand'].unique():

        group.loc[group['start'].isin(ex_1start_more_ends['start']), 'splicing_event'] ='alt_5_prime'
        group.loc[group['end'].isin(ex_1end_more_starts['end']), 'splicing_event'] ='alt_3_prime'
        if len(composite_condition1)>0:
            group.loc[group['start'].isin(composite_condition1['start']), 'splicing_event'] ='composite_alt_5_skipped'
        if len(composite_condition2)>0:
            group.loc[group['end'].isin(composite_condition2['end']), 'splicing_event'] ='composite_alt_3_skipped'
        group.loc[group.apply(lambda row: (row['start'], row['end']) in skipped_tuples and pd.isna(row['splicing_event']), axis=1),
        'splicing_event'] = 'skipped_exon'
        group.loc[group.apply(lambda row: ((row['start'], row['end']) in constitutive_tuples) and ((row['start'], row['end']) in constitutive_alltrs_tuples), axis=1),
        'splicing_event'] = 'constitutive'
        group.loc[group.apply(lambda row: ((row['start'], row['end']) in constitutive_tuples) and not ((row['start'], row['end']) in constitutive_alltrs_tuples), axis=1),
        'splicing_event'] = 'constitutive_totranscript'

    elif '-' in group['strand'].unique():
        group.loc[group['start'].isin(ex_1start_more_ends['start']), 'splicing_event'] ='alt_3_prime(-)'
        if len(composite_condition1)>0:
            group.loc[group['start'].isin(composite_condition1['start']), 'splicing_event'] ='composite_alt_3_skipped(-)'
        group.loc[group['end'].isin(ex_1end_more_starts['end']), 'splicing_event'] ='alt_5_prime(-)'
        if len(composite_condition2)>0:
            group.loc[group['end'].isin(composite_condition2['end']), 'splicing_event'] ='composite_alt_5_skipped(-)'
        group.loc[group.apply(lambda row: (row['start'], row['end']) in skipped_tuples and pd.isna(row['splicing_event']), axis=1),
        'splicing_event'] = 'skipped_exon(-)'
        group.loc[group.apply(lambda row: ((row['start'], row['end']) in constitutive_tuples) and ((row['start'], row['end']) in constitutive_alltrs_tuples), axis=1),
        'splicing_event'] = 'constitutive(-)'
        group.loc[group.apply(lambda row: ((row['start'], row['end']) in constitutive_tuples) and not ((row['start'], row['end']) in constitutive_alltrs_tuples), axis=1),
        'splicing_event'] = 'constitutive_totranscript(-)'
    # This snippet is repeated so that some exons, which get labelled wrongly as any other events, even though they show those characteristics, to keep it simple, this will override and will always label such exons as first exons over any other thing.
    if '+' in group['strand'].unique():
        group.loc[smallest_exon_rows.index, 'splicing_event'] = 'first_exon'
        group.loc[greatest_exon_rows.index, 'splicing_event'] = 'last_exon'

    elif '-' in group['strand'].unique():
        group.loc[smallest_exon_rows.index, 'splicing_event'] = 'last_exon'
        group.loc[greatest_exon_rows.index, 'splicing_event'] = 'first_exon'

    gene_model = gene_model.merge(
        group[['start', 'end', 'splicing_event']],
        on=['start', 'end'],
        how='left'
    )
    gene_model['splicing_event'] = gene_model['splicing_event_y'].combine_first(gene_model['splicing_event_x'])
    gene_model.drop(['splicing_event_x', 'splicing_event_y'], axis=1, inplace=True)

    return gene_model

def intron_filter(group):
    exon_condition = group[
    group['splicing_event'].astype(str).str.contains('constitutive') |
    group['splicing_event'].astype(str).str.contains('skipped_exon') | 
    group['splicing_event'].astype(str).str.contains('first') |
    group['splicing_event'].astype(str).str.contains('last')]
    exon_range = set(exon_condition[['start','end']].itertuples(index=False, name=None))
    first_exons = group[group['splicing_event'].astype(str).str.contains('first')]
    last_exons  = group[group['splicing_event'].astype(str).str.contains('last')]
    alt_3_exons = group[group['splicing_event'].astype(str).str.contains('alt_3')]
    alt_5_exons = group[group['splicing_event'].astype(str).str.contains('alt_5')]
    skipped_exons = group[group['splicing_event'].astype(str).str.contains('skipped_exon')]
    skipped_range = set(skipped_exons[['start', 'end']].itertuples(index=False, name=None))
    constitutive = group[group['splicing_event'].astype(str).str.contains('constitutive')]
    
    if '+' in group['strand'].unique():
        limitcoord_alt3 = alt_3_exons.groupby('end').apply(
            lambda group: group.loc[group['start'].idxmin()]
        ).reset_index(drop=True)

        limitcoord_alt5 = alt_5_exons.groupby('start').apply(
            lambda group: group.loc[group['end'].idxmax()]
        ).reset_index(drop=True)
        
    elif '-'  in group['strand'].unique():
        limitcoord_alt5 = alt_5_exons.groupby('end').apply(
            lambda group: group.loc[group['start'].idxmin()]
        ).reset_index(drop=True)

        limitcoord_alt3 = alt_3_exons.groupby('start').apply(
            lambda group: group.loc[group['end'].idxmax()]
        ).reset_index(drop=True)
    
    group['correct_coords'] = np.nan
    introns = group[group['group']=='intron'].copy()
    
    for index, row in introns.iterrows():
        start, end, strand = row['start'], row['end'], row['strand']
        tr_id = row['transcript_id']
        is_in_exon_range =  any(
            (exon_start <= end and exon_end >= start)
            for exon_start, exon_end in exon_range)
        
        if strand == '+':
            matches_first = ((start in (first_exons['end'] + 1).values) or (end in (first_exons['start'] - 1).values)) & (tr_id in first_exons['transcript_id'])
            
            matches_last = ((start in (last_exons['end'] + 1).values) or (end in (last_exons['start'] - 1).values)) & (tr_id in last_exons['transcript_id'])
            
            matches_alt3 = not limitcoord_alt3.empty and (
                end in (limitcoord_alt3['start'] - 1).values)
            
            matches_alt5 = not limitcoord_alt5.empty and (
                 start in (limitcoord_alt5['end'] + 1).values)

            matches_constitutive = not constitutive.empty and (
                start in (constitutive['end'] + 1).values or end in (constitutive['start'] - 1).values)
            
            matches_skipped = not skipped_exons.empty and (
                start in (skipped_exons['end'] + 1).values or end in (skipped_exons['start'] - 1).values)
        
        elif (strand == '-'):
            matches_first = ((start in (first_exons['end'] + 1).values) or (end in (first_exons['start'] - 1).values)) & (tr_id in first_exons['transcript_id'])

            matches_last = ((start in (last_exons['end'] + 1).values) or (end in (last_exons['start'] - 1).values)) & (tr_id in last_exons['transcript_id'])
            
            matches_alt3 = not limitcoord_alt3.empty and (
                start in (limitcoord_alt3['end'] + 1).values)
            
            matches_alt5 = not limitcoord_alt5.empty and (
                end in (limitcoord_alt5['start'] - 1).values)

            matches_constitutive = not constitutive.empty and (
                    start in (constitutive['end'] + 1).values or end in (constitutive['start'] - 1).values)

            matches_skipped = not skipped_exons.empty and (
                    start in (skipped_exons['end'] + 1).values or end in (skipped_exons['start'] - 1).values)
            
            
        if (matches_alt3 or matches_alt5 or matches_constitutive or matches_skipped or matches_first or matches_last) and not is_in_exon_range:
            group.at[index, 'correct_coords'] = 'True'
        else:
            group.at[index, 'correct_coords'] = None

    for index, row in group.iterrows():
        type = row['group']
        if type == 'intron' and pd.isna(row['correct_coords']):
            group.drop(index, inplace=True)

    return group


def intron_cat(group):
    introns = group[group['group']=='intron']
    first_int = group[group['splicing_event'].str.contains('first', na = False)][['start','end','strand','transcript_id']]
    last_int = group[group['splicing_event'].str.contains('last', na = False)][['start','end','strand','transcript_id']]

    for index, row in introns.iterrows():
        if  (row['strand']=='+'):
            if (row['start'] in first_int['end'].values + 1) and (row['transcript_id'] in first_int['transcript_id'].values):
                current_val = group.at[index, 'splicing_event']
                group.at[index, 'splicing_event'] = f"{current_val}, first_intron" if pd.notna(current_val) else "first_intron"
            elif (row['end'] in last_int['start'].values - 1) and (row['transcript_id'] in last_int['transcript_id'].values):
                current_val = group.at[index, 'splicing_event']
                group.at[index, 'splicing_event'] = f"{current_val}, last_intron" if pd.notna(current_val) else "last_intron"
        elif  (row['strand']=='-'):         
            if (row['end'] in first_int['start'].values - 1) and (row['transcript_id'] in first_int['transcript_id'].values):
                current_val = group.at[index, 'splicing_event']
                group.at[index, 'splicing_event'] = f"{current_val}, first_intron(-)" if pd.notna(current_val) else "first_intron(-)"
            elif (row['start'] in last_int['end'].values + 1) and (row['transcript_id'] in last_int['transcript_id'].values):
                current_val = group.at[index, 'splicing_event']
                group.at[index, 'splicing_event'] = f"{current_val}, last_intron(-)" if pd.notna(current_val) else "last_intron(-)"

    constex_1 = group[group['splicing_event'].str.contains('constitutive',na=False)][['start','end','strand']]

    # Labels introns flanking constitutive exons on both strands
    for index, row in introns.iterrows():
        if (row['start'] in constex_1['end'].values + 1):
            if  (row['strand']=='+'):
                group.at[index, 'splicing_event'] = 'constitutive_downstream'    
            elif (row['strand']=='-'):                                                                          
                group.at[index, 'splicing_event'] = 'constitutive_upstream(-)'               
                                                                                                                 
        elif (row['end'] in constex_1['start'].values - 1):
            if  (row['strand']=='+'):
                group.at[index, 'splicing_event'] = 'constitutive_upstream'    
            elif (row['strand']=='-'):                                                               
                group.at[index, 'splicing_event'] = 'constitutive_downstream(-)'

    constex = group[group['splicing_event'].str.contains('constitutive_totranscript',na=False)][['start','end','strand']]

    # Labels introns flanking constitutive exons grouped by first exon coords, on both strands
    for index, row in introns.iterrows():
        if (row['start'] in constex['end'].values + 1):
            if  (row['strand']=='+'):
                group.at[index, 'splicing_event'] = 'constitutive_totranscript_downstream'    
            elif (row['strand']=='-'):                                                                          
                group.at[index, 'splicing_event'] = 'constitutive_totranscript_upstream(-)'               
                                                                                                                 
        elif (row['end'] in constex['start'].values - 1):
            if  (row['strand']=='+'):
                group.at[index, 'splicing_event'] = 'constitutive_totranscript_upstream'    
            elif (row['strand']=='-'):                                                               
                group.at[index, 'splicing_event'] = 'constitutive_totranscript_downstream(-)'

    # Labels introns flanking skipped exons on both strands   
    skipex = group[group['splicing_event'].str.contains('skipped_exon',na=False)][['start','end','strand']]
    for index, row in introns.iterrows():
        if (row['start'] in skipex['end'].values + 1):
            if (row['strand']=='+'):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', skipped_downstream'
                else:
                    group.at[index, 'splicing_event'] = 'skipped_downstream'
            elif (row['strand']=='-'):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', skipped_upstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'skipped_upstream(-)'
            
        elif (row['end'] in skipex['start'].values - 1):
            if (row['strand']=='+'):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', skipped_upstream'
                else:
                    group.at[index, 'splicing_event'] = 'skipped_upstream'
            elif (row['strand']=='-'):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', skipped_downstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'skipped_downstream(-)'
                    
        #Labels alternative 5'ss exon flanking introns on both strands:
    alt_5 = group[group['splicing_event'].str.contains('alt_5_prime',na=False)][['gene_id','group','start','end','strand','transcript_id']]
    alt_5_ex = alt_5[alt_5['group']=='exon']
    if '+' in alt_5['strand'].unique():
        alt5_ends = alt_5_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()
        max_end_indices = alt5_ends.groupby(['gene_id', 'start'])['end'].idxmax() #Now I am only selecting the start and end coordinates which are #containing a greater end position i.e. the second 5'ss Because this is where the intron truly begins.
        alt5_grt_end = alt5_ends.loc[max_end_indices]
    elif '-' in alt_5['strand'].unique():
        alt_5_starts = alt_5_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()
        min_start_indices = alt_5_starts.groupby(['gene_id', 'end'])['start'].idxmin()
        neg_alt_small_starts = alt_5_starts.loc[min_start_indices]

    for index, row in introns.iterrows():
        if (row['strand']=='+') and not alt_5.empty:
            if (row['start'] in alt5_grt_end['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_5_downstream'
                else:
                    group.at[index, 'splicing_event'] = 'alt_5_downstream'
            elif (row['end'] in alt5_grt_end['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_5_upstream'
                else:
                    group.at[index, 'splicing_event'] = 'alt_5_upstream'
                        
        elif (row['strand']=='-') and not alt_5.empty:
            if (row['start'] in neg_alt_small_starts['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_5_upstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'alt_5_upstream(-)'
            elif (row['end'] in neg_alt_small_starts['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_5_downstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'alt_5_downstream(-)'
                    
    #Labelling composite alternative 5'ss skipped exon flanking introns 
    composite_alt_5 = group[group['splicing_event'].str.contains('composite_alt_5_skipped',na=False)][['gene_id','group','start','end','strand','transcript_id']]
    composite_alt_5_ex = composite_alt_5[composite_alt_5['group']=='exon']
    if '+' in composite_alt_5['strand'].unique():
        comp_alt5_ends = composite_alt_5_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()
        comp_max_end_indices = comp_alt5_ends.groupby(['gene_id', 'start'])['end'].idxmax() #Now I am only selecting the start and end coordinates which are #containing a greater end position i.e. the second 5'ss Because this is where the intron truly begins.
        comp_alt5_grt_end = comp_alt5_ends.loc[comp_max_end_indices]
    elif '-' in composite_alt_5['strand'].unique():
        composite_alt_5_starts = composite_alt_5_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()
        comp_min_start_indices = composite_alt_5_starts.groupby(['gene_id', 'end'])['start'].idxmin()
        comp_neg_alt_small_starts = composite_alt_5_starts.loc[comp_min_start_indices]

    for index, row in introns.iterrows():
        if (row['strand']=='+') and not composite_alt_5.empty:
            if (row['start'] in comp_alt5_grt_end['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_5_downstream'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_5_downstream'
            elif (row['end'] in comp_alt5_grt_end['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_5_upstream'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_5_upstream'

        elif (row['strand']=='-') and not composite_alt_5.empty:
            if (row['start'] in comp_neg_alt_small_starts['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_5_upstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_5_upstream(-)'
            elif (row['end'] in comp_neg_alt_small_starts['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_5_downstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_5_downstream(-)'

    #Labels alternative 5'ss exon flanking introns on both strands:
    alt_3 = group[group['splicing_event'].str.contains('alt_3_prime',na=False)][['gene_id','group','start','end','strand','transcript_id']]
    alt_3_ex = alt_3[alt_3['group']=='exon']
    if '+' in alt_3['strand'].unique():
        alt3_starts = alt_3_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()
        min_start_indices = alt3_starts.groupby(['gene_id', 'end'])['start'].idxmin()
        alt3_small_starts = alt3_starts.loc[min_start_indices]
    elif '-' in alt_3['strand'].unique():
        N_alt3_ends = alt_3_ex.groupby(['gene_id', 'start', 'end'])['transcript_id'].nunique().reset_index()
        max_end_indices = N_alt3_ends.groupby(['gene_id', 'start'])['end'].idxmax()
        N_alt3_grt_end = N_alt3_ends.loc[max_end_indices]


    for index, row in introns.iterrows():
        if (row['strand']=='+') and not alt_3.empty:
            if (row['start'] in alt3_small_starts['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_3_downstream'
                else:
                    group.at[index, 'splicing_event'] = 'alt_3_downstream'
            elif (row['end'] in alt3_small_starts['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_3_upstream'
                else:
                    group.at[index, 'splicing_event'] = 'alt_3_upstream'

        elif (row['strand']=='-') and not alt_3.empty:
            if (row['start'] in N_alt3_grt_end['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_3_upstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'alt_3_upstream(-)'
            elif (row['end'] in N_alt3_grt_end['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', alt_3_downstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'alt_3_downstream(-)'

    #Labelling composite alternative 3'ss skipped exon flanking introns 
    composite_alt_3 = group[group['splicing_event'].str.contains('composite_alt_3_skipped',na=False)][['gene_id','group','start','end','strand','transcript_id']]
    composite_alt_3_ex = composite_alt_3[composite_alt_3['group']=='exon']
    if '+' in composite_alt_3['strand'].unique():
        comp_alt3_starts = composite_alt_3_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()
        comp_min_start_indices = comp_alt3_starts.groupby(['gene_id', 'end'])['start'].idxmin()
        comp_alt3_small_starts = comp_alt3_starts.loc[comp_min_start_indices]
    elif '-' in composite_alt_3['strand'].unique():
        comp_N_alt3_ends = composite_alt_3_ex.groupby(['gene_id', 'start', 'end'])['transcript_id'].nunique().reset_index()
        comp_max_end_indices = comp_N_alt3_ends.groupby(['gene_id', 'start'])['end'].idxmax()
        comp_N_alt3_grt_end = comp_N_alt3_ends.loc[comp_max_end_indices]


    for index, row in introns.iterrows():
        if (row['strand']=='+') and not composite_alt_3.empty:
            if (row['start'] in comp_alt3_small_starts['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_3_downstream'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_3_downstream'
            elif (row['end'] in comp_alt3_small_starts['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_3_upstream'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_3_upstream'

        elif (row['strand']=='-') and not composite_alt_3.empty:
            if (row['start'] in comp_N_alt3_grt_end['end'].values + 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_3_upstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_3_upstream(-)'
            elif (row['end'] in comp_N_alt3_grt_end['start'].values - 1):
                current_val = group.at[index, 'splicing_event']
                if pd.notna(current_val):
                    group.at[index, 'splicing_event'] = current_val + ', composite_alt_3_downstream(-)'
                else:
                    group.at[index, 'splicing_event'] = 'composite_alt_3_downstream(-)'
                   
    return group
        

keywords = []

import pandas as pd

r2 = pd.read_csv(gff_file, sep='\t', names=['Chr', 'source', 'group', 'start', 'end', '_', 'strand', '-', 'info'])
r2.dropna(axis=0, inplace=True)
r2_Ex_Int = r2[(r2['group'] == 'exon') | (r2['group'] == 'intron')].copy()

# r2_Ex_Int[['ID','Parent','gene_id', 'gene_symbol','transcript_id','transcript_symbol']] = r2_Ex_Int[
# 'info'].str.split(';', expand=True) #- Drosophila gff has these entries in metadata column

print('\n \\ The original gff dataframe. Parsing out the metadata column - info. \\\n')
print(r2_Ex_Int.head(10))

with open(keywords_file, 'r') as file:
    for w in file:
        keywords.append(str(w.rstrip()))

# Converting metadata column identifiers (keywords) into columns and extracting relevant information.
pattern = '|'.join(r'(?P<{}>{}=[^;]+)'.format(k, k) for k in keywords)
extracted_info = r2_Ex_Int['info'].str.extractall(pattern)
grouped_info = extracted_info.fillna('').astype(str).groupby(level=0).agg(lambda x: ';'.join(filter(None, x)))
grouped_info.columns = [f'{col}' for col in grouped_info.columns]
r2_Ex_Int = r2_Ex_Int.join(grouped_info, how='left')

print('\n \\ Here\'s a glimpse of the modified gff dataframe. \\ \n')
print(r2_Ex_Int.head(10))
r2_Ex_Int.to_csv('gff_parsed_meta.csv', sep='\t')


r2_Ex_Int.sort_values(by=['gene_id', 'start'], ascending=True, inplace=True)
r2_Ex_Int.drop(labels='info', axis=1, inplace=True)
r2_Ex_Int.loc[:, 'sort_order'] = r2_Ex_Int.apply(lambda x: x['start'] if x['strand'] == '+' else -x['start'], axis=1)
r2_Ex_Int = r2_Ex_Int.sort_values(by=['sort_order'])
r2_Ex_Int = r2_Ex_Int.drop(columns=['sort_order'])

r2_Ex_Int['order_of_appearance'] = r2_Ex_Int.groupby(['group', 'transcript_id']).cumcount() + 1

r2_Ex_Int_1 = r2_Ex_Int.groupby(['gene_id']).apply(lambda x: retained_introns(x))
r2_Ex_Int_1.reset_index(drop=True, inplace=True)
r2_Exons_annot_composite = r2_Ex_Int_1.groupby(['gene_id']).apply(lambda x: exon_cat(x))
r2_Exons_annot_composite.reset_index(drop=True, inplace=True)

retained_introns = len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'retained intron'])
constitutive = len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'constitutive']) + len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'constitutive(-)'])
skipped = len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'skipped_exon']) + len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'skipped_exon(-)'])
alt_5 = len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'alt_5_prime']) + len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'alt_5_prime(-)'])
alt_3 = len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'alt_3_prime']) + len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'alt_3_prime(-)'])
composite = len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'composite_alt_3_skipped']) + len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'composite_alt_3_skipped(-)']) + len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'composite_alt_5_skipped']) + len(r2_Exons_annot_composite[r2_Exons_annot_composite['splicing_event'] == 'composite_alt_5_skipped(-)'])

print(f'Total number of constitutive exons: {constitutive}')
print(f'Total number of skipped exons: {skipped}')
print(f'Total number of alternative 5\'ss exons: {alt_5}')
print(f'Total number of alternative 3\'ss exons: {alt_3}')
print(f'Total number of composite alternative 5\'ss or alternative 3\'ss containing and skipped exons: {composite}')
print(f'Total number of retained introns: {retained_introns}')
df_final = r2_Exons_annot_composite.groupby(['gene_id']).apply(lambda x: intron_filter(x))
df_final.reset_index(drop=True, inplace=True)
all_annot_df = df_final.groupby(['gene_id']).apply(lambda x: intron_cat(x))

df_final.to_csv('./Exonslabelled_intronsnot.csv', sep='\t')

print(all_annot_df.head(10))

all_annot_df.reset_index(drop=True, inplace=True)

all_annot_df.sort_values(by=['gene_id', 'transcript_id', 'order_of_appearance'], ascending=True, inplace=True)

############################################################# The functions below, scaled_coordinates () and scale_group () were written for cases where we only have sequence fasta files of orthologues and reference sps
# For example, the analysis I did on C.elegans switch like splicing events. I had relative coordinates that were inferred from the fasta files. I didn't have absolute coordinates ######## If this isn't the case, the lines  
# where these functions are being called can be commented out. ##### The lines are only to represent everything in terms of sclaed or relative coordinates.################################################################# 
if coordinate_type == 'relative':
    def scale_coordinates(row, start, end):# removed strand from args, since same function works for both strands.
        max_length = start - end
        qstart = 0
        scaled_rstart = qstart + (row['start'] - start)
        scaled_rend = qstart + (row['end'] - start)
        return pd.Series([scaled_rstart, scaled_rend], index=['scaled_rstart', 'scaled_rend'])

        # Group by 'gene' and apply scaling function
    def scale_group(group):
        # For each group (gene), calculate min and max start and end, works for both strands since start and ends are always 
        # smaller start and greater end, only the order in negative strand is reversed, i.e., exon with smallest start is the 
        # last exon. Draw it out to avoid confusion.
        strand = group['strand'].iloc[0]
        min_rstart = group['start'].min()
        max_rend = group['end'].max()
        group[['scaled_rstart', 'scaled_rend']] = group.apply(scale_coordinates, axis=1, args=(min_rstart, max_rend))
        
        group = group.rename(columns={'start': 'absolute_start', 'end': 'absolute_end'})
        group = group.rename(columns={'scaled_rstart': 'start', 'scaled_rend': 'end'})
        return group

    # Assuming 'Parent' identifies different genes
    scaled_df = all_annot_df.groupby('Parent').apply(scale_group).reset_index(drop=True)
    # Assuming 'Parent' identifies different genes
    all_annot_Introns = scaled_df[scaled_df['group'] == 'intron']

elif coordinate_type =='absolute':
    all_annot_Introns = all_annot_df[all_annot_df['group'] == 'intron']

all_annot_Introns = all_annot_Introns.copy()
#Unknown_Intron_indices = all_annot_Introns[all_annot_Introns['splicing_event_y_y'] == 'intron'].index
#all_annot_Introns.drop(Unknown_Intron_indices, inplace=True)
all_annot_Introns['ID'] = all_annot_Introns['Chr'].astype(str) + '_' + all_annot_Introns['start'].astype(str) + '_' + \
                          all_annot_Introns['end'].astype(str)
all_annot_Introns.groupby(['ID', 'gene_id', 'splicing_event'])['splicing_event'].nunique()
Intron_df = all_annot_Introns

# 1. Calculating size of introns
all_annot_Introns['size'] = abs(all_annot_Introns['start'] - all_annot_Introns['end'])
print('The median size of introns is : ')
print(all_annot_Introns['size'].median())

median_size = all_annot_Introns['size'].median()


def round_to_nearest_multiple_of_10(number):
    return round(number / 10) * 10

if median == 0:
    median = round_to_nearest_multiple_of_10(median_size)

print('The introns bigger than the median size are being resized to - %d bp. \n' % median)


# 2.  Eliminating extremely small introns
small_introns = all_annot_Introns[all_annot_Introns['size'] < 10].index
all_annot_Introns = all_annot_Introns.drop(small_introns)

# 3. Creating a separate dataframe for introns smaller than median bp in length
Intron_size_median = all_annot_Introns[all_annot_Introns['size'] < median]

# 4. Processing introns greater than median bp in length
Intron_size_great_median = all_annot_Introns[all_annot_Introns['size'] >= median]
Intron_size_great_median['End_medianbp_apart'] = Intron_size_great_median['start'] + median
Intron_size_great_median['Start_medianbp_apart'] = Intron_size_great_median['end'] - median

# Creating new ID columns for starting median bp and ending median bp of introns
Intron_size_great_median['ID_Starting_median'] = Intron_size_great_median['Chr'].astype(str) + '_' + \
                                                 Intron_size_great_median[ 
                                                     'strand'] + '_' + Intron_size_great_median['start'].astype(
    str) + '_' + Intron_size_great_median[
                                                     'End_medianbp_apart'].astype(str)

Intron_size_great_median['ID_Ending_median'] = Intron_size_great_median['Chr'].astype(str) + '_' + \
                                               Intron_size_great_median[
                                                   'strand'] + '_' + Intron_size_great_median[
                                                   'Start_medianbp_apart'].astype(str) + '_' + Intron_size_great_median[
                                                   'end'].astype(str)

merged_intron_df = pd.concat([Intron_size_great_median, Intron_size_median], ignore_index=True)
merged_intron_df.drop_duplicates(inplace=True)
merged_intron_df.to_csv(intron_annot_file, sep='\t')

# Extracting the coordinates and arranging them based on their strands. ##
# NOTE: THESE NEED TO BE CONVERTED TO BED FILE ###
# CAN BE CONVERTED TO BED FILE EASILY WITH EXCEL AND R - GTF-to-bed.r###
intron_beginningmedian = Intron_size_great_median[
    ['Parent','Chr','group', 'start', 'End_medianbp_apart', 'strand', 'ID_Starting_median']]
intron_endingmedian = Intron_size_great_median[['Parent','Chr','group', 'Start_medianbp_apart', 'end', 'strand', 'ID_Ending_median']]
small_intron_coordinates = Intron_size_median[['Parent','Chr','group', 'start', 'end', 'strand', 'ID']]
intron_beginningmedian.rename(columns={'End_medianbp_apart': 'end', 'ID_Starting_median': 'ID'}, inplace=True)
intron_endingmedian.rename(columns={'Start_medianbp_apart': 'start', 'ID_Ending_median': 'ID'}, inplace=True)
big_introns = pd.concat([intron_beginningmedian, intron_endingmedian], ignore_index=True)
all_intron_coordinates = pd.concat([small_intron_coordinates, big_introns], ignore_index=True)
all_intron_coordinates.drop_duplicates(inplace=True)
all_intron_coordinates.to_csv(intron_cords_file, sep='\t')


if coordinate_type == 'relative':
    all_exon_coordinates = scaled_df[scaled_df['group']=='exon'][['Parent','Chr','group','start', 'end','strand', 'ID']]
    all_coordinates = pd.concat([all_intron_coordinates, all_exon_coordinates], ignore_index=True)
    all_coordinates.sort_values(by=['Parent', 'start'], ascending=[True, True], inplace=True)
    
    all_coordinates.to_csv('./CelsWS15_Exon_IntronFragment_coordinates.tsv', sep='\t', index=False)
    scaled_df.to_csv(exon_intron_annot_file, sep='\t')

elif coordinate_type == 'absolute':
    all_exon_coordinates = all_annot_df[all_annot_df['group']=='exon'][['Parent','Chr','group','start', 'end','strand', 'ID']]
    all_coordinates = pd.concat([all_intron_coordinates, all_exon_coordinates], ignore_index=True)
    all_coordinates.sort_values(by=['Parent', 'start'], ascending=[True, True], inplace=True)

    all_coordinates.to_csv('./Exon_IntronFragment_coordinates.tsv', sep='\t', index=False)
    all_annot_df.to_csv(exon_intron_annot_file, sep='\t')
# The end - this will give an annotated intron coordinates file which can be further analyzed in excel.###

#%%
