# Inferring exon and intron annotations from genome feature file (.gff)

## gff_exon-intron_annotations.py for processing .gff -

1. DOWNLOADING GTF/GFF file -
The .gff fle can be downloaded from RefSeq or other specific databases.
These usually contain the coordinates of a gene and all other entities within it such as exons, UTRs, CDS etc.
These files usually lack introns hence, intron coordinates need to be inferred from the exon coordinates.

2. USING AGAT software TO INFER INTRON COORDINATES-
It can be downloaded and installed from this source -  https://github.com/NBISweden/AGAT
This can be further added to the path variable and following script can be run
agat_sp_add_introns.pl --gff [_.gff] -o [_updated.gff]

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

1. metadata_Introns_annotated.tsv -> This file only contains actual intron-coordinates and inferred metadata.Depending on the options chosen [absolute/relative] [size], columns will be added to the file.
   1.a. If [absolute] is chosen - start and end coordinates will be present in start and end columns
   1.b. If [relative] is chosen - start and end coordinates will be present in absolute start and absolute_end columns while the relative coordinates are 	calculated assuming each gene starts at 0.
   1.c. [size] - Will be 'chosen size - median or supplied size' and less bp
   1.d. IDs of intron fragments with a certain size will also be present in columns
   
3. metadata_Introns_Exons_annotated.tsv -> This file contains both introns and exons related metadata.

5. intron_coordinate.tsv -> This file contains correct inferred coordinates of introns. This file can be further converted to bed file using. crdnts_to_bed.R 

GENERATING BED FILES FROM COORDINATES FILE -

The R script - crdnts_to_bed.R takes the coordinates.tsv file generated previously as input.
It then outputs two bed files. One is for + strand sequences and another is for - strand sequences.

The bed file is generated using convert2bed function but the bed files are not tab delimited. Hence, these positive_stand.bed (+) and
negative_strand.bed (-) files can be manually edited in a vi editor or other editor of choice to change the spaces to tabs.

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
