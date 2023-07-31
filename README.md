# ortholog_seq_extract
The CDS, intron and intergenic sequences from each one-to-one ortholog pair are extracted and aligned using blastn. The output  is the alignment-length-weighted similarity.
## running example
input files needed: genome and annotation files from the two genomes. a one-to-one orthologs table <br />

python3 ortholog_seq_extract.py genome_1.fa,genome_2.fa genome_1.gff,genome_2.gff one-to-one_orthologs.tsv result.tsv
