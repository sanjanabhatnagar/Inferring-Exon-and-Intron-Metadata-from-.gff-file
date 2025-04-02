# Inferring Exon and Intron Annotations from Genome Feature Files (.gff)

## Overview
This pipeline extracts exon and intron metadata from .gff files using only coordinate information. It classifies exons into the following splicing categories:

- **Constitutively spliced exons** (always included in mRNA)
- **Skipped exons** (e.g., casette exons and mutually exclusive exons)
- **Alternative 5’ splice site (5’ss) exons** (exons with two donor sites)
- **Alternative 3’ splice site (3’ss) exons** (exons with two acceptor sites)
- **Retained introns** (introns incorporated into mRNA)
- **Composite exons** (exhibiting multiple alternative splicing types)

Additionally, it annotates the first and last exons per transcript for a given gene model.

## Inferring Intron Coordinates
The .gff file does not include intron information. To address this:
1. **AGAT** software infers inter-exonic region coordinates.
2. **gff_exon-intron_annotations.py** filters actual introns and adds metadata using flanking exon annotations.
3. The **crdnts_to_bed.R** script converts extracted intron coordinates to BED format.

---

## Pipeline Steps
### 1. Download the .gff File
.gff files can be obtained from RefSeq or other genomic databases. They typically include gene structure elements (exons, UTRs, CDS, etc.) but lack intron coordinates.

### 2. Infer Intron Coordinates Using AGAT
[AGAT](https://github.com/NBISweden/AGAT) can be installed and added to the system PATH. Use the following command:
```bash
agat_sp_add_introns.pl --gff input.gff -o updated.gff
```
This step generates an updated .gff file with inter-exonic regions.

### 3. Annotate Exons and Introns
Run the Python script to process the updated .gff file:
```bash
python gff_exon-intron_annotations.py updated.gff keywords.txt metadata_Introns_annotated.tsv metadata_Introns_Exons_annotated.tsv intron_coordinates.tsv [absolute/relative] [size]
```
#### Parameters:
- **updated.gff**: Output from AGAT containing inferred intron regions.
- **keywords.txt**: List of metadata fields to extract (e.g., `ID`, `gene_id`, `transcript_id`). Since, information column in .gff for different organisms contains different data.
- **absolute/relative**: Determines coordinate type. `absolute` outputs standard genomic coordinates, while `relative` scales them from 0 to gene length.
- **size**: Defines intron fragment length (default = median intron size, rounded to the nearest 10 bp).

#### Intron Splitting Strategy:
- Introns larger than the median size are split into two parts:
  - **First** fragment (start of the intron)
  - **Last** fragment (end of the intron)
  
---

## Output Files
### 1. `metadata_Introns_annotated.tsv`
- Contains actual intron coordinates and metadata.
- If `relative` is chosen, absolute and relative start/end coordinates are included.
- IDs for intron fragments are generated.

### 2. `metadata_Introns_Exons_annotated.tsv`
- Includes metadata for both exons and introns.

### 3. `intron_coordinates.tsv`
- Contains filtered intron coordinates, which can be converted to BED format.

---

## Converting to BED Format
Use `crdnts_to_bed.R` to convert `intron_coordinates.tsv` to BED format:
```r
Rscript crdnts_to_bed.R intron_coordinates.tsv
```
This generates two BED files:
- **Positive strand (`positive_strand.bed`)**
- **Negative strand (`negative_strand.bed`)**

### Adding Unique IDs
Use `awk` to append IDs:
```bash
awk '{printf("%s_+_%d_%d\n",$1,$2,$3); }' positive_strand.bed > positive_ID.bed
paste positive_strand.bed positive_ID.bed > final_positive.bed

awk '{printf("%s_-_%d_%d\n",$1,$2,$3); }' negative_strand.bed > negative_ID.bed
paste negative_strand.bed negative_ID.bed > final_negative.bed
```
After editing to replace spaces with tabs, the output should resemble:
```
NC_026501.1     3822    3894    NC_026501.1_+_3822_3894
NC_026501.1     4861    4961    NC_026501.1_+_4861_4961
...
```
---

## Summary
This pipeline extracts exon and intron annotations from .gff files by inferring intron positions, filtering actual introns, and generating BED files for downstream analysis. It supports metadata extraction, intron classification, and size-based intron fragmentation.

For any issues, refer to AGAT documentation or modify the scripts accordingly.

## References 

Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  
(Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717

