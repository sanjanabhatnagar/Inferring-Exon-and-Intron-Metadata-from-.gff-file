# Takes approximately 11 minutes to run. ##
import os
import sys

gtf_file = sys.argv[1]  # This is the input gtf
keywords_file = sys.argv[2]  # This is the column labels we want the gtf metadata - "info" column to be parsed into.
# Reason - Every gtf's info column consists different identifier (types/categories of information) and we need the
# identifiers to become columns in final files.
intron_annot_file = sys.argv[3]  # This is the gtf with all exons and introns annotated and metadata
intron_cords_file = sys.argv[4]  # This is just the coordinates of introns ready to be converted to bed format

# Depending on the length you want your sequences to be resized to, you can either choose to enter a length here
# alternatively if you want to use the median length for resizing you can comment this out and use the code below on
# line 399 : median = round_to_nearest_multiple_of_10(median_size)
# To keep it simple, I have used the variable name median and it depends on user to use the actual median size or
#  specify a length they'd like the introns to be resized to.
median = sys.argv[5]
# CAN BE CONVERTED TO BED FILE EASILY WITH EXCEL AND R - GTF-to-bed.r ###

keywords = []

import pandas as pd

# The big gtf file with all groups - transcripts, exons, UTRs etc.
#
r2 = pd.read_csv(gtf_file, sep='\t', names=['Chr', 'source', 'group', 'start', 'end', '_', 'strand', '-', 'info'])
r2.dropna(axis=0, inplace=True)
r2_Ex_Int = r2[(r2['group'] == 'exon') | (r2['group'] == 'intron')].copy()

# r2_Ex_Int[['ID','Parent','gene_id', 'gene_symbol','transcript_id','transcript_symbol']] = r2_Ex_Int[
# 'info'].str.split(';', expand=True) #- Drosophila gff has these entries in metadata column

print('\n \\ The original gff dataframe. Parsing out the metadata column - info. \\\n')
print(r2_Ex_Int.head(10))

# The following code reads the info column (the column holding metadata)
with open(keywords_file, 'r') as file:
    for w in file:
        keywords.append(str(w.rstrip()))

# Converting metadata column identifiers (keywords) into columns and extracting relevant information.
pattern = '|'.join(r'(?P<{}>{}=[^;]+)'.format(k, k) for k in keywords)

# Extracts information from pattern and converts into columns.
extracted_info = r2_Ex_Int['info'].str.extractall(pattern)

grouped_info = extracted_info.fillna('').astype(str).groupby(level=0).agg(lambda x: ';'.join(filter(None, x)))
# Rename columns of grouped information DataFrame
grouped_info.columns = [f'{col}' for col in grouped_info.columns]
# Join the grouped columns with the original DataFrame
r2_Ex_Int = r2_Ex_Int.join(grouped_info, how='left')

print('\n \\ Here\'s a glimpse of the modified gff dataframe. \\ \n')
print(r2_Ex_Int.head(10))
r2_Ex_Int.to_csv('gff_parsed_meta.csv', sep='\t')

# Basic sorting of dataframe based on certain columns
r2_Ex_Int.sort_values(by=['gene_id', 'start'], ascending=True, inplace=True)
r2_Ex_Int.drop(labels='info', axis=1, inplace=True)

r2_Ex_Int.loc[:, 'sort_order'] = r2_Ex_Int.apply(lambda x: x['start'] if x['strand'] == '+' else -x['start'], axis=1)
r2_Ex_Int = r2_Ex_Int.sort_values(by=['sort_order'])
r2_Ex_Int = r2_Ex_Int.drop(columns=['sort_order'])

# Adding order of appearance column for different exons and introns
r2_Ex_Int['order_of_appearance'] = r2_Ex_Int.groupby(['group', 'transcript_id']).cumcount() + 1

# # 1. Classifying exons based on splicing events (4 classes - constitutive, alternative 5'ss, alternative 3'ss and
# skipped exons) Main snippet which classifies exons
r2_Exons = r2_Ex_Int[r2_Ex_Int['group'] == 'exon']
r2_Exons['splicing_event'] = ''

r2_Exons_P = r2_Ex_Int[r2_Ex_Int['strand'] == '+']
r2_Exons_N = r2_Ex_Int[r2_Ex_Int['strand'] == '-']

# Class 1 - Identify alternative 5'ss containing exons
P_alt_5_prime_exons = r2_Exons_P.groupby(['gene_id', 'start']).apply(
    lambda group: group if group['end'].nunique() > 1 else pd.DataFrame()).droplevel(0)
transcript_ids_1 = set(P_alt_5_prime_exons.index.get_level_values(0).tolist())
r2_Exons.loc[r2_Exons['start'].isin(transcript_ids_1), 'splicing_event'] = 'alt_5_prime'

N_alt_5_prime_exons = r2_Exons_N.groupby(['gene_id', 'end']).apply(
    lambda group: group if group['start'].nunique() > 1 else pd.DataFrame()).droplevel(0)
transcript_ids_4 = set(N_alt_5_prime_exons.index.get_level_values(0).tolist())
r2_Exons.loc[r2_Exons['end'].isin(transcript_ids_4), 'splicing_event'] = 'alt_5_prime'

# Class 2 - Identify alternative 3'ss containing exons
P_alt_3_prime_exons = r2_Exons_P.groupby(['gene_id', 'end']).apply(
    lambda group: group if group['start'].nunique() > 1 else pd.DataFrame()).droplevel(0)
transcript_ids_2 = set(P_alt_3_prime_exons.index.get_level_values(0).tolist())
r2_Exons.loc[r2_Exons['end'].isin(transcript_ids_2), 'splicing_event'] = 'alt_3_prime'

N_alt_3_prime_exons = r2_Exons_N.groupby(['gene_id', 'start']).apply(
    lambda group: group if group['end'].nunique() > 1 else pd.DataFrame()).droplevel(0)
transcript_ids_3 = set(N_alt_3_prime_exons.index.get_level_values(0).tolist())
r2_Exons.loc[r2_Exons['start'].isin(transcript_ids_3), 'splicing_event'] = 'alt_3_prime'

# Class 3 - Identify Constitutive exons 3.1 Grouping by gene_id, start, and end, then counting the unique
# transcript_ids so that I know how many exons belong to a particular gene
exon_counts = r2_Exons.groupby(['gene_id', 'start', 'end'])['transcript_id'].count()

# 3.2 The following line tells me how many transcripts belong to a particular gene
total_unique_transcripts = r2_Exons.groupby('gene_id')['transcript_id'].nunique().reset_index()

# 3.3 Merging on common columns to align indices
merged_counts = pd.merge(exon_counts.reset_index(), total_unique_transcripts, on='gene_id', how='left',
                         suffixes=('_count', '_unique'))

# 3.4 Selecting rows where the counts match - The logic I am using is if the (unique) number of an exon == (unique)
# number of transcripts Which means that the exon is present in all transcripts for a gene
constitutive_exons = merged_counts[merged_counts['transcript_id_count'] == merged_counts['transcript_id_unique']]

# 3.5 Labelling them as constitutive exons
r2_Exons.loc[r2_Exons.set_index(['gene_id', 'start', 'end']).index.isin(
    constitutive_exons.set_index(['gene_id', 'start', 'end']).index), 'splicing_event'] = 'Constitutive'

# Class 4 - Identify skipped exons

r2_Exons_rem = r2_Exons[r2_Exons['splicing_event'] == '']

rem_exon_counts = r2_Exons_rem.groupby(['gene_id', 'start', 'end'])['transcript_id'].count()

# 4.1 Checking if the count of transcripts is equal to the count of unique transcript_ids
rem_total_unique_transcripts = r2_Exons_rem.groupby('gene_id')['transcript_id'].nunique().reset_index()

rem_merged_counts = pd.merge(rem_exon_counts.reset_index(), rem_total_unique_transcripts, on='gene_id', how='left',
                             suffixes=('_count', '_unique'))

# 4.2 Selecting rows where the counts of (unique) exons is less than the (unique) number of transcripts which means
# this exon is not present in all transcripts of a gene
skipped_exons = rem_merged_counts[rem_merged_counts['transcript_id_count'] < rem_merged_counts['transcript_id_unique']]

# 4.3 Labelling them as skipped exons
r2_Exons.loc[r2_Exons.set_index(['gene_id', 'start', 'end']).index.isin(
    skipped_exons.set_index(['gene_id', 'start', 'end']).index), 'splicing_event'] = 'Skipped_Exon'

# Getting a new dataframe that has this dataframe and introns
df_merged = pd.merge(r2_Ex_Int, r2_Exons, how='left')
df_merged['splicing_event'] = df_merged['splicing_event'].fillna('intron')

# 2. Inferring itron types based on exon events
# 2.1. Labelling Constitutive exon flanking introns (+-)
# Runtime 400-600 ms
# 1. Segmenting the dataset - Only selecting +ve strand entries
pos = df_merged[(df_merged['strand'] == '+')]

# 2. Further selecting constitutive exons in a separate dataframe
constex = pos[pos['splicing_event'].str.contains('Constitutive')]
# 3. Now getting all the start and end coordinates for constitutive exons/exons present in all transcripts of a gene.
constex_df = constex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()

# 4. Finding flanking introns.
downstream_intron_condition = pos[(pos['group'] == 'intron') & pos['start'].isin(constex_df['end'] + 1)].index
upstream_intron_condition = pos[(pos['group'] == 'intron') & pos['end'].isin(constex_df['start'] - 1)].index

# 5. Updating the 'splicing_event' column based on conditions
pos.loc[downstream_intron_condition, 'splicing_event'] = 'Downstr_intron-Constitutive'
pos.loc[upstream_intron_condition, 'splicing_event'] = 'Upstr_intron-Constitutive'

# 1. Segmenting the dataset - Only selecting -ve strand entries
neg = df_merged[(df_merged['strand'] == '-')]

# 2. Further selecting constitutive exons in a separate dataframe
constex_neg = neg[neg['splicing_event'].str.contains('Constitutive')]

# 3. Now getting all the start and end coordinates for constitutive exons/exons present in all transcripts of a gene.
constex_neg_df = constex_neg.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()

# 4. Finding flanking introns.
upstream_intron_neg_condition = neg[(neg['group'] == 'intron') & neg['start'].isin(constex_neg_df['end'] - 1)].index
downstream_intron_neg_condition = neg[(neg['group'] == 'intron') & neg['end'].isin(constex_neg_df['start'] + 1)].index

# Update the 'splicing_event' column based on conditions
neg.loc[upstream_intron_neg_condition, 'splicing_event'] = 'Upstr_intron_neg-Constitutive'
neg.loc[downstream_intron_neg_condition, 'splicing_event'] = 'Downstr_intron_neg-Constitutive'

# 2.2. Labelling Alternative 5'ss exon flanking introns (+-)
# Runtime 108 ms
# 1. Segmenting the dataset - Only selecting +ve strand alt_5_prime events both introns and exons

# 2. Then extracting exons so that I can work on their coordinates to figure out the right end coordinate in this case
# right coordinate out of two end coordinates as this will be used to infer the correct intron start coordinate


pos_alt5_ex = pos[pos['splicing_event'].str.contains('alt_5_prime')]
# 2.1 Now getting all the start and end coordinates for this alt_5 ss containing exons.
# The following line gives the number of transcripts containing an exon with certain start and end coordinates.
alt5_ends = pos_alt5_ex.groupby(['gene_id', 'start', 'end'])['transcript_id'].nunique().reset_index()

# 2.2 Now I am only selecting the start and end coordinates which are containing a greater end position i.e. the
# second 5'ss Because this is where the intron truly begins. Using this I will find the true intron boundary and
# correct coordinates for downstream intron.
max_end_indices = alt5_ends.groupby(['gene_id', 'start'])['end'].idxmax()
alt5_grt_end = alt5_ends.loc[max_end_indices]

# 3. Now I iterate over the rows of the dataframe containing both exons and introns found by AGAT software to find the true introns.
# AGAT software, illogically, infers introns considering both 5ss of preceding exon which means some inferred introns, specially the ones
# inferred from the first 5'ss will be containing some portion of exon, up until the second 5ss which wasn't used in the transcript.

UI_alt5_pos = pos[(pos['group'] == 'intron') & pos['end'].isin(alt5_grt_end['start'] - 1)].index
DI_alt5_pos = pos[(pos['group'] == 'intron') & pos['start'].isin(alt5_grt_end['end'] + 1)].index

# Update the 'splicing_event' column based on conditions
pos.loc[UI_alt5_pos, 'splicing_event'] = 'Upstr_Intron-alt_5_prime'
pos.loc[DI_alt5_pos, 'splicing_event'] = 'Downstr_Intron-alt_5_prime'

# 1. Segmenting the dataset - Only selecting -ve strand neg_alt5__prime events both introns and exons


# 2. Then extracting exons so that I can work on their coordinates to figure out the right start coordinate in this case
# right coordinate out of two start coordinates as this will be used to infer the correct intron end coordinate
neg_alt5_ex = neg[neg['splicing_event'].str.contains('alt_5_prime')]

##2.1 Now getting all the start and end coordinates for this neg_alt_ ss containing exons.
# The following line gives the number of transcripts containing an exon with certain start and end coordinates.
neg_alt_starts = neg_alt5_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()

##2.2 Now I am only selecting the start and end coordinates which are containing a smaller start negition i.e. the first 'ss
## Because this is where the intron truly ends.
## Using this I will find the true intron boundary and correct coordinates for upstream intron.
min_start_indices = neg_alt_starts.groupby(['gene_id', 'end'])['start'].idxmin()
neg_alt_small_starts = neg_alt_starts.loc[min_start_indices]

UI_alt5_neg = neg[(neg['group'] == 'intron') & neg['start'].isin(neg_alt_small_starts['end'] + 1)].index
DI_alt5_neg = neg[(neg['group'] == 'intron') & neg['end'].isin(neg_alt_small_starts['start'] - 1)].index

# Update the 'splicing_event' column based on conditions
neg.loc[UI_alt5_neg, 'splicing_event'] = 'Upstr_Intron_neg-alt_5_prime'
neg.loc[DI_alt5_neg, 'splicing_event'] = 'Downstr_Intron_neg-alt_5_prime'

# . Now I iterate over the rows of the dataframe containing both exons and introns found by AGAT software to find the true introns.
# AGAT software, illogically, nonsensically infers introns considering both ss of suceeding exon which means some inferred introns, specially the ones
# inferred from the first 'ss will be containing some portion of exon, up until the second ss which wasn't used in the transcript.

## 2.3. Labelling Alternative 3'ss exon flanking introns (+-)

pos2 = pos
# 1. Segmenting the dataset - Only selecting +ve strand alt_3_prime events both introns and exon
pos_alt3_ex = pos2[pos2['splicing_event'].str.contains('alt_3_prime')]

# 2. Then extracting exons so that I can work on their coordinates to figure out the right start coordinate in this case

##2.1 Now getting all the start and end coordinates for this alt_3 ss containing exons.
# The following line gives the number of transcripts containing an exon with certain start and end coordinates.
alt3_starts = pos_alt3_ex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()

##2.2 Now I am only selecting the start and end coordinates which are containing a smaller start position i.e. the first 3'ss
## Because this is where the intron truly ends.
## Using this I will find the true intron boundary and correct coordinates for upstream intron.
min_start_indices = alt3_starts.groupby(['gene_id', 'end'])['start'].idxmin()
alt3_small_starts = alt3_starts.loc[min_start_indices]

# 3. Now I iterate over the rows of the dataframe containing both exons and introns found by AGAT software to find the true introns.
# AGAT software, illogically, nonsensically infers introns considering both 3ss of suceeding exon which means some inferred introns, specially the ones
# inferred from the first 3'ss will be containing some portion of exon, up until the second 3ss which wasn't used in the transcript.

UI_alt3_pos = pos2[(pos2['group'] == 'intron') & pos2['end'].isin(alt3_small_starts['start'] - 1)].index
DI_alt3_pos = pos2[(pos2['group'] == 'intron') & pos2['start'].isin(alt3_small_starts['end'] + 1)].index

pos2.loc[UI_alt3_pos, 'splicing_event'] = 'Upstr_Intron-alt_3_prime'
pos2.loc[DI_alt3_pos, 'splicing_event'] = 'Downstr_Intron-alt_3_prime'

neg2 = neg
# 1. Segmenting the dataset - Only selecting +ve strand alt_3_prime events both introns and exon
N_alt3_ex = neg2[neg2['splicing_event'].str.contains('alt_3_prime')]
N_alt3_ends = N_alt3_ex.groupby(['gene_id', 'start', 'end'])['transcript_id'].nunique().reset_index()

##2.2 Now I am only selecting the start and end coordinates which are containing a greater end position i.e. the second 5'ss
## Because this is where the intron truly begins.
## Using this I will find the true intron boundary and correct coordinates for downstream intron.
max_end_indices = N_alt3_ends.groupby(['gene_id', 'start'])['end'].idxmax()
N_alt3_grt_end = N_alt3_ends.loc[max_end_indices]

DI_alt3_neg = neg2[(neg2['group'] == 'intron') & neg2['end'].isin(N_alt3_grt_end['start'] - 1)].index
UI_alt3_neg = neg2[(neg2['group'] == 'intron') & neg2['start'].isin(N_alt3_grt_end['end'] + 1)].index

neg2.loc[UI_alt3_neg, 'splicing_event'] = 'Upstr_Intron_neg-alt_3_prime'
neg2.loc[DI_alt3_neg, 'splicing_event'] = 'Downstr_Intron_neg-alt_3_prime'

## Now merging both alt_5_prime +-  and alt_3_prime +- dfs


col_names = ['Chr', 'group', 'start', 'end', 'strand', 'order_of_appearance'] + keywords
pos_df = pd.merge(pos, pos2, how='outer',
                  on=col_names)
# Group by the relevant columns and aggregate the 'splicing_event' values
pos_df = pos_df.groupby(col_names).agg({
    'splicing_event_y': lambda x: ', '.join(x.dropna())
}).reset_index()

neg_df = pd.merge(neg, neg2, how='outer', on=col_names)
# Group by the relevant columns and aggregate the 'splicing_event' values
neg_df = neg_df.groupby(col_names).agg({
    'splicing_event_y': lambda x: ', '.join(x.dropna())
}).reset_index()

#  Final dataframe with all annotations, none overwritten!!
final_df = pd.concat([pos_df, neg_df], ignore_index=True)

#  2.4. Labelling Skipped exon flanking introns (+-)
# 1. Segmenting the dataset - Only selecting +ve strand skipped exon events both introns and exons
fdf_1 = final_df
# 1. Segmenting the dataset - Only selecting +ve strand skipped exons.
skipex = fdf_1[(fdf_1['splicing_event_y'].str.contains('Skipped_Exon')) & (fdf_1['strand'] == '+')]

# 2.1 Now getting all the start and end coordinates for skipped exons.
# The following line gives the number of transcripts containing an exon with certain start and end coordinates.
skipex_df = skipex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()

DI_SE_pos = fdf_1[(fdf_1['group'] == 'intron') & fdf_1['start'].isin(skipex_df['end'] + 1)].index
UI_SE_pos = fdf_1[(fdf_1['group'] == 'intron') & fdf_1['end'].isin(skipex_df['start'] - 1)].index

fdf_1.loc[DI_SE_pos, 'splicing_event_y'] = 'Downstr_intron-skipped_exon'
fdf_1.loc[UI_SE_pos, 'splicing_event_y'] = 'Upstr_intron-skipped_exon'

# 1. Segmenting the dataset - Only selecting -ve strand skipped exons.
skipex_neg = fdf_1[(fdf_1['splicing_event_y'].str.contains('Skipped_Exon')) & (fdf_1['strand'] == '-')]

# 2.1 Now getting all the start and end coordinates for skipped exons.
# The following line gives the number of transcripts containing an exon with certain start and end coordinates.
skipex_neg_df = skipex.groupby(['gene_id', 'end', 'start'])['transcript_id'].nunique().reset_index()

UI_SE_neg = fdf_1[(fdf_1['group'] == 'intron') & fdf_1['start'].isin(skipex_neg_df['end'] + 1)].index
DI_SE_neg = fdf_1[(fdf_1['group'] == 'intron') & fdf_1['end'].isin(skipex_neg_df['start'] - 1)].index

fdf_1.loc[DI_SE_neg, 'splicing_event_y'] = 'Downstr_intron_neg-skipped_exon'
fdf_1.loc[UI_SE_neg, 'splicing_event_y'] = 'Upstr_intron_neg-skipped_exon'

all_annot_df = pd.merge(fdf_1, final_df, how='outer', on=col_names)
# Group by the relevant columns and aggregate the 'splicing_event' values
all_annot_df = all_annot_df.groupby(col_names).agg({
    'splicing_event_y_y': lambda x: ', '.join(x.dropna())
}).reset_index()


# 3. Filtering out any introns that still may contain exons (wrongly inferred by AGAT) + sanity check Function to
# check if an intron's start and end coordinates contain any exon's start and end coordinates. This is the error
# introduced  by AGAT software when it infers introns from exon coordinates

def contains_exon(intron_starts, intron_ends, exon_starts, exon_ends, strand):
    exon_starts = exon_starts.values.reshape((-1, 1))
    exon_ends = exon_ends.values.reshape((-1, 1))
    intron_starts = intron_starts.values
    intron_ends = intron_ends.values
    if strand == '-':
        return ((intron_starts <= exon_ends) & (intron_ends >= exon_starts)).any(axis=0)
    elif strand == '+':
        return ((intron_starts >= exon_starts) & (intron_ends <= exon_ends)).any(axis=0)


# Filter out introns based on exon containment
filtered_introns = []
introns_to_remove_indices = []
for gene_id, group in all_annot_df.groupby('gene_id'):
    strand = group['strand'].iloc[0]
    exons = group[group['group'] == 'exon']
    introns = group[group['group'] == 'intron']
    exon_starts = exons['start']
    exon_ends = exons['end']
    intron_starts = introns['start']
    intron_ends = introns['end']
    introns_to_discard = contains_exon(intron_starts, intron_ends, exon_starts, exon_ends, strand)
    filtered_introns.append(introns[introns_to_discard])
    if introns_to_discard.any():
        introns_to_remove_indices.extend(introns.index[introns_to_discard])

# Concatenate the filtered introns back into a DataFrame
filtered_introns_df = pd.concat(filtered_introns)

all_annot_df.drop(introns_to_remove_indices, inplace=True)
all_annot_df.sort_values(by=['gene_id', 'transcript_id', 'order_of_appearance'], ascending=True, inplace=True)

all_annot_Introns = all_annot_df[all_annot_df['group'] == 'intron']

all_annot_Introns = all_annot_Introns.copy()
Unknown_Intron_indices = all_annot_Introns[all_annot_Introns['splicing_event_y_y'] == 'intron'].index
all_annot_Introns.drop(Unknown_Intron_indices, inplace=True)
all_annot_Introns['ID'] = all_annot_Introns['Chr'].astype(str) + '_' + all_annot_Introns['start'].astype(str) + '_' + \
                          all_annot_Introns['end'].astype(str)
all_annot_Introns.groupby(['ID', 'gene_id', 'splicing_event_y_y'])['splicing_event_y_y'].nunique()

# The order of exon classes is in a way that I label constitutive exon flanking introns first so that the labels
# inferred based on alternatively spliced exons take precedence. For example there might be certain downstream
# constitutive introns which are upstream of skipped exon so we want the label to be Upstr_intron-skipped_exon
# instead. 4. Resizing introns - For Introns greater than median bp in size (taking first and last median bp) and
# keeping the smaller introns in a separate dataframe

Intron_df = all_annot_Introns

# 1. Calculating size of introns
all_annot_Introns['size'] = abs(all_annot_Introns['start'] - all_annot_Introns['end'])
print('The median size of introns is : ')
print(all_annot_Introns['size'].median())

median_size = all_annot_Introns['size'].median()


def round_to_nearest_multiple_of_10(number):
    return round(number / 10) * 10


#median = round_to_nearest_multiple_of_10(median_size)
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
    ['Chr', 'start', 'End_medianbp_apart', 'strand', 'ID_Starting_median']]
intron_endingmedian = Intron_size_great_median[['Chr', 'Start_medianbp_apart', 'end', 'strand', 'ID_Ending_median']]
small_intron_coordinates = Intron_size_median[['Chr', 'start', 'end', 'strand', 'ID']]
intron_beginningmedian.rename(columns={'End_medianbp_apart': 'end', 'ID_Starting_median': 'ID'}, inplace=True)
intron_endingmedian.rename(columns={'Start_medianbp_apart': 'start', 'ID_Ending_median': 'ID'}, inplace=True)
big_introns = pd.concat([intron_beginningmedian, intron_endingmedian], ignore_index=True)
all_intron_coordinates = pd.concat([small_intron_coordinates, big_introns], ignore_index=True)
all_intron_coordinates.drop_duplicates(inplace=True)
all_intron_coordinates.to_csv(intron_cords_file, sep='\t')

# The end - this will give an annotated intron coordinates file which can be further analyzed in excel.###

#%%
