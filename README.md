# Inferring exon and intron annotations from genome feature file (.gff)

## GTF_ExInt_Annotator.py for Processing gff -

1. DOWNLOADING GTF/GFF file -
The GFF fle can be downloaded from RefSeq or other specific databases.
These usually contain the coordinates of a gene and all other entities within it such as exons, UTRs, CDS etc.
These files usually lack introns hence, intron coordinates need to be inferred from the exon coordinates.

2. USING AGAT software TO INFER INTRON COORDINATES-
It can be downloaded and installed from this source -  https://github.com/NBISweden/AGAT
This can be further added to the path variable and following script can be run
agat_sp_add_introns.pl --gff [_.gff] -o [_introns.gff]

3. CLASSIFYING INTRONS -
The script gff_exon-intron_annotations.py specifically annotates all the exons and further infers intron coordinates and annotations based on the exons they flank. The usage is -

python gff_exon-intron_annotations.py [_updated.gff] [keywords.txt] [metadata_Introns_annotated.tsv] [metadata_Introns_Exons_annotated.tsv][intron_coordinates.tsv] [absolute/relative] [size]

[_updated.gff] - This is the file generated previously using AGAT software which adds the coordinates of regions between exons across different transcripts in a gene model to the existing gff. Not all regions depict actual introns!! Hence, the gff_exon-intron_annotations.py has a function that filters actual intron coordinates from this updated file.

[keywords.txt] - This depends on the metadata column of a gff. Every gff is different in terms of the information it contains in the 'info' column.
The identifiers which the user wants to extract data for can be supplied in the text file. For example -

```
ID
Parent
gene_id
gene_symbol
transcript_id
transcript_symbol

```
Only the values associated with these keywords in the metadata column will be extracted and will be stored into additional
columns.
The script uses 'gene_id', start and end coordinates and transcript_ids for classification so these keywords are must.

In order to maintain a uniformity in sizes of the introns being generated i.e. if they are to be subjected to downstream
analysis, this script calculates the median intron size and rounds this up to a closest multiple of 10 and takes all the
introns smaller than this and the introns bigger than the median size are broken down into two parts. However, if the user
Wants to specify a certain length they can enter the size in the end as indicated in the command above.

The starting region of Introns until a new end coordinate which is median-rounded-to-closest-multiple-of-10 (or user-specified length) bp apart.
The ending region of Introns but a new start coordinate which is median-rounded-to-closest-multiple-of-10 (or user-specified length) bp upstream of
end coordinate.
For example if median intron size is 96, introns bigger than this will be cropped to 100 in way that starting 100 bp of
this intron (from donor 5'splice site end) and ending 100 bp of this intron (from acceptor 3'splice site end) will be taken.

[absolute/relative] - This is the type of coordinates for the output intron coordinates file. In all cases, the output will be 'absolute.' The 'relative' option is used only when the output file requires gene coordinates scaled from 0 to the full length of the gene.

3.1. RESULTING FILES -

metadata_Introns_annotated.tsv -> This file only contains introns and related metadata

	Chr	source	group	absolute_start	absolute_end	_	strand	-	ID	Parent	exon_id	exon_number	gene_biotype	gene_id	gene_source	transcript_biotype	transcript_id	transcript_source	order_of_appearance	splicing_event	correct_coords	start	end	size	End_medianbp_apart	Start_medianbp_apart	ID_Starting_median	ID_Ending_median
0	II	WormBase	intron	1912	2505	.	+	.	II_64_657	Parent=2L52.1a.1	exon_id=2L52.1a.1.e7	exon_number=7	gene_biotype=protein_coding	gene_id=WBGene00007063	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2L52.1a.1	transcript_source=WormBase	1	first_intron, skipped_upstream	TRUE	64	657	593	134	587	II_+_64_134	II_+_587_657
1	II	WormBase	intron	3037	3405	.	+	.	II_1189_1557	Parent=2L52.1a.1	exon_id=2L52.1a.1.e7	exon_number=7	gene_biotype=protein_coding	gene_id=WBGene00007063	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2L52.1a.1	transcript_source=WormBase	4	skipped_downstream	TRUE	1189	1557	368	1259	1487	II_+_1189_1259	II_+_1487_1557
2	II	WormBase	intron	3553	3801	.	+	.	II_1705_1953	Parent=2L52.1a.1	exon_id=2L52.1a.1.e7	exon_number=7	gene_biotype=protein_coding	gene_id=WBGene00007063	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2L52.1a.1	transcript_source=WormBase	5	first_intron, constitutive_upstream, skipped_downstream	TRUE	1705	1953	248	1775	1883	II_+_1705_1775	II_+_1883_1953
3	II	WormBase	intron	3985	4200	.	+	.	II_2137_2352	Parent=2L52.1a.1	exon_id=2L52.1a.1.e7	exon_number=7	gene_biotype=protein_coding	gene_id=WBGene00007063	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2L52.1a.1	transcript_source=WormBase	6	last_intron, constitutive_downstream	TRUE	2137	2352	215	2207	2282	II_+_2137_2207	II_+_2282_2352
4	II	WormBase	intron	3553	3801	.	+	.	II_17_265	Parent=2L52.1b.1	exon_id=2L52.1b.1.e3	exon_number=3	gene_biotype=protein_coding	gene_id=WBGene00007063	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2L52.1b.1	transcript_source=WormBase	1	first_intron, constitutive_upstream, skipped_downstream	TRUE	17	265	248	87	195	II_+_17_87	II_+_195_265
5	II	WormBase	intron	3985	4200	.	+	.	II_449_664	Parent=2L52.1b.1	exon_id=2L52.1b.1.e3	exon_number=3	gene_biotype=protein_coding	gene_id=WBGene00007063	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2L52.1b.1	transcript_source=WormBase	2	last_intron, constitutive_downstream	TRUE	449	664	215	519	594	II_+_449_519	II_+_594_664
6	II	WormBase	intron	15268555	15269458	.	+	.	II_343_1246	Parent=2RSSE.1a.1	exon_id=2RSSE.1a.1.e5	exon_number=5	gene_biotype=protein_coding	gene_id=WBGene00007064	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2RSSE.1a.1	transcript_source=WormBase	1	first_intron, constitutive_upstream	TRUE	343	1246	903	413	1176	II_+_343_413	II_+_1176_1246
7	II	WormBase	intron	15270032	15270795	.	+	.	II_1820_2583	Parent=2RSSE.1a.1	exon_id=2RSSE.1a.1.e5	exon_number=5	gene_biotype=protein_coding	gene_id=WBGene00007064	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2RSSE.1a.1	transcript_source=WormBase	3	last_intron, constitutive_downstream, skipped_upstream	TRUE	1820	2583	763	1890	2513	II_+_1820_1890	II_+_2513_2583
8	II	WormBase	intron	15268555	15269458	.	+	.	II_345_1248	Parent=2RSSE.1b.1	exon_id=2RSSE.1b.1.e8	exon_number=8	gene_biotype=protein_coding	gene_id=WBGene00007064	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2RSSE.1b.1	transcript_source=WormBase	1	first_intron, constitutive_upstream	TRUE	345	1248	903	415	1178	II_+_345_415	II_+_1178_1248
9	II	WormBase	intron	15270032	15270795	.	+	.	II_1822_2585	Parent=2RSSE.1b.1	exon_id=2RSSE.1b.1.e8	exon_number=8	gene_biotype=protein_coding	gene_id=WBGene00007064	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2RSSE.1b.1	transcript_source=WormBase	3	last_intron, constitutive_downstream, skipped_upstream	TRUE	1822	2585	763	1892	2515	II_+_1822_1892	II_+_2515_2585
10	II	WormBase	intron	15274517	15275028	.	+	.	II_6307_6818	Parent=2RSSE.1b.1	exon_id=2RSSE.1b.1.e8	exon_number=8	gene_biotype=protein_coding	gene_id=WBGene00007064	gene_source=WormBase	transcript_biotype=protein_coding	transcript_id=2RSSE.1b.1	transcript_source=WormBase	5	skipped_downstream	TRUE	6307	6818	511	6377	6748	II_+_6307_6377	II_+_6748_6818![image](https://github.com/user-attachments/assets/2cd889fb-7b44-400e-8369-35adffc738e6)


metadata_Introns_Exons_annotated.tsv -> This file contains both introns and exons related metadata, if interested in looking at coding region annotations.
intron_coordinate.tsv -> This file contains correct inferred coordinates of introns belonging to different categories. 
This file can be further converted to bed file using 
GTF-to-bed.R 

GENERATING BED FILES FROM COORDINATES FILE -

The R script - GTF-to-bed.R takes the coordinates.tsv file generated previously as input.
It then outputs two bed files. One is for + strand sequences and another is for - strand sequences.

The bed file is generated using convert2bed function but the bed files are not tab delimited. Hence, these pv.bed (+) and
nv.bed (-) files can be manually edited in a vi editor to change the spaces to tabs.

Then in order to add IDs which are lost after running the R script, new ID column can be made using awk command.

```
awk '{printf("%s_+_%d_%d\n"),$1,$2,$3,$4; }' pv.bed > pv_2.bed # I added + in ID for positive strand seqs
paste pv.bed pv_2.bed > pv_ID.bed

#for negative strand seqs -
awk '{printf("%s_-_%d_%d\n"),$1,$2,$3,$4; }' nv.bed > nv_2.bed
paste nv.bed nv_2.bed > nv_ID.bed
```
Then the file can be viewed and first row has to be omitted. The first row contains column labels which can result in errors later on.

It should be looking like this

```
NC_026501.1     3822    3894    NC_026501.1_+_3822_3894
NC_026501.1     4861    4961    NC_026501.1_+_4861_4961
NC_026501.1     5013    5113    NC_026501.1_+_5013_5113
NC_026501.1     12189   12244   NC_026501.1_+_12189_12244
NC_026501.1     16016   16064   NC_026501.1_+_16016_16064
NC_026501.1     16296   16396   NC_026501.1_+_16296_16396
NC_026501.1     16390   16490   NC_026501.1_+_16390_16490
NC_026501.1     18363   18437   NC_026501.1_+_18363_18437
NC_026501.1     34876   34976   NC_026501.1_+_34876_34976

```
