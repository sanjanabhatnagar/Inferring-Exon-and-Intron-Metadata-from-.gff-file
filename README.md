# Understanding RNA grammar : Functional grouping of non-coding regions in high-dimensional RNA-binding protein motif space

## GTF_ExInt_Annotator.py for Processing gff -

1. DOWNLOADING GTF/GFF file -
The GFF fle can be downloaded from RefSeq or other specific databases.
These usually contain the coordinates of a gene and all other entities within it such as exons, UTRs, CDS etc.
These files usually lack introns hence, intron coordinates need to be inferred from the exon coordinates.

2. USING AGAT software TO INFER INTRON COORDINATES-
It's pretty easy to install and can be found here -  https://github.com/NBISweden/AGAT
This can be further added to the path variable and following script can be run
agat_sp_add_introns.pl --gff [_.gff] -o [_introns.gff]

3. CLASSIFYING INTRONS -
The script GTF_ExInt_Annotator.py specifically annotates all the inferred introns based on the exons they flank. The usage is -

python GTF_ExInt_Annotator.py [_introns.gff] [keywords.txt] [metadata_containing_annotated.tsv] [intron_coordinates.tsv] [size]

[_introns.gff] - This is the file generated previously using AGAT software which adds the introns to the existing gff.
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
Wants to specify a certain length they can comment out line 399 : median = round_to_nearest_multiple_of_10(median_size) in the script.
And specify the length with the python command
python GTF_ExInt_Annotator.py [_introns.gff] [keywords.txt] [metadata_containing_annotated.tsv] [intron_coordinates.tsv] [size]

The starting region of Introns until a new end coordinate which is median-rounded-to-closest-multiple-of-10 (or user-specified length) bp apart.
The ending region of Introns but a new start coordinate which is median-rounded-to-closest-multiple-of-10 (or user-specified length) bp upstream of
end coordinate.
For example if median intron size is 96, introns bigger than this will be cropped to 100 in way that starting 100 bp of
this intron and ending 100 bp of this intron will be taken.

3.1. RESULTING FILES -

metadata_containing_annotated.tsv -> This file only contains introns and related metadata
gff_parsed_meta.tsv -> This file contains exons and introns and related metadata
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
