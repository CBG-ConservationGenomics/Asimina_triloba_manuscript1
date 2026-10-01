#!/bin/bash

# ============================================================
# Download genome, protein, and annotation data
#
# Genomes:
#   genome_sequences/
#
# Original/unfiltered proteins:
#   unfiltered_protein_sequences/
#
# Annotations:
#   annotations/
#
# Raw downloaded files:
#   downloads/
# ============================================================


# Create directories if they do not already exist
mkdir -p downloads
mkdir -p genome_sequences
mkdir -p annotations
mkdir -p unfiltered_protein_sequences


echo "========================================"
echo "Downloading genome and protein data"
echo "========================================"
echo


# ============================================================
# Annona cherimola – UMA
# ============================================================

echo "=== Annona cherimola ==="

# Genome
wget -O downloads/Annona_cherimola.genome.fa.gz \
  "https://ihsmsubtropicals.uma.es/downloads/Annona%20cherimola/Sequences/Anche102_genome.fasta.gz"

gunzip -c downloads/Annona_cherimola.genome.fa.gz \
  > genome_sequences/Annona_cherimola.fa


# Protein – primary transcripts only
wget -O downloads/Annona_cherimola.protein.faa.gz \
  "https://ihsmsubtropicals.uma.es/downloads/Annona%20cherimola/Sequences/anche102_proteins_primaryTranscriptOnly_annot.fasta.gz"

gunzip -c downloads/Annona_cherimola.protein.faa.gz \
  > unfiltered_protein_sequences/Annona_cherimola.fa


echo "Annona cherimola complete!"
echo


# ============================================================
# Annona montana – CNCB
# ============================================================

echo "=== Annona montana ==="

# Genome
wget -O downloads/Annona_montana.genome.fa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Annona_montana_Am_v1.0_GWHDQZG00000000/GWHDQZG00000000.genome.fasta.gz"

gunzip -c downloads/Annona_montana.genome.fa.gz \
  > genome_sequences/Annona_montana.fa


# Protein
wget -O downloads/Annona_montana.protein.faa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Annona_montana_Am_v1.0_GWHDQZG00000000/GWHDQZG00000000.Protein.faa.gz"

gunzip -c downloads/Annona_montana.protein.faa.gz \
  > unfiltered_protein_sequences/Annona_montana.fa


# GFF annotation
wget -O downloads/Annona_montana.gff.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Annona_montana_Am_v1.0_GWHDQZG00000000/GWHDQZG00000000.gff.gz"

gunzip -c downloads/Annona_montana.gff.gz \
  > annotations/Annona_montana.gff


echo "Annona montana complete!"
echo


# ============================================================
# Lindera megaphylla – CNCB
# ============================================================

echo "=== Lindera megaphylla ==="

# Genome
wget -O downloads/Lindera_megaphylla.genome.fa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Lindera_megaphylla_LMv1_GWHBKHA00000000/GWHBKHA00000000.genome.fasta.gz"

gunzip -c downloads/Lindera_megaphylla.genome.fa.gz \
  > genome_sequences/Lindera_megaphylla.fa


# Protein
wget -O downloads/Lindera_megaphylla.protein.faa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Lindera_megaphylla_LMv1_GWHBKHA00000000/GWHBKHA00000000.Protein.faa.gz"

gunzip -c downloads/Lindera_megaphylla.protein.faa.gz \
  > unfiltered_protein_sequences/Lindera_megaphylla.fa


# GFF annotation
wget -O downloads/Lindera_megaphylla.gff.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Lindera_megaphylla_LMv1_GWHBKHA00000000/GWHBKHA00000000.gff.gz"

gunzip -c downloads/Lindera_megaphylla.gff.gz \
  > annotations/Lindera_megaphylla.gff


echo "Lindera megaphylla complete!"
echo


# ============================================================
# Magnolia kwangsiensis – CNCB
# ============================================================

echo "=== Magnolia kwangsiensis ==="

# Genome
wget -O downloads/Magnolia_kwangsiensis.genome.fa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Magnolia_kwangsiensis_Mkfd_GWHGEUP00000000.1/GWHGEUP00000000.1.genome.fasta.gz"

gunzip -c downloads/Magnolia_kwangsiensis.genome.fa.gz \
  > genome_sequences/Magnolia_kwangsiensis.fa


# Protein
wget -O downloads/Magnolia_kwangsiensis.protein.faa.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Magnolia_kwangsiensis_Mkfd_GWHGEUP00000000.1/GWHGEUP00000000.1.Protein.faa.gz"

gunzip -c downloads/Magnolia_kwangsiensis.protein.faa.gz \
  > unfiltered_protein_sequences/Magnolia_kwangsiensis.fa


# GFF annotation
wget -O downloads/Magnolia_kwangsiensis.gff.gz \
  "https://download.cncb.ac.cn/gwh/Plants/Magnolia_kwangsiensis_Mkfd_GWHGEUP00000000.1/GWHGEUP00000000.1.gff.gz"

gunzip -c downloads/Magnolia_kwangsiensis.gff.gz \
  > annotations/Magnolia_kwangsiensis.gff


echo "Magnolia kwangsiensis complete!"
echo


# ============================================================
# Persea americana – downloaded manually
# ============================================================

# Files downloaded manually:
#   pame.pep
#   pame.fa.gz
#   pame.gff.gz
#
# These should ultimately be organized as:
#
#   genome_sequences/Persea_americana.fa
#   unfiltered_protein_sequences/Persea_americana.fa
#   annotations/Persea_americana.gff


echo "=== Persea americana downloaded manually ==="
echo


# ============================================================
# Asimina triloba – local files
# ============================================================

# Asimina triloba files are available locally and therefore
# are not downloaded by this script.
#
# They should ultimately be organized as:
#
#   genome_sequences/Asimina_triloba.fa
#   unfiltered_protein_sequences/Asimina_triloba.fa
#   annotations/Asimina_triloba.gff


echo "=== Asimina triloba is local (not downloaded here) ==="
echo


# ============================================================
# Species excluded from final analysis
# ============================================================

# Annona muricata:
# Protein sequences are scaffold-based and the corresponding
# genome/GFF needed to identify chromosome-anchored proteins
# are not available.
#
# Cinnamomum micranthum:
# Assembly consists of scaffold-level sequences and a reliable
# chromosome/pseudomolecule protein subset could not be defined.


# ============================================================
# Done
# ============================================================

echo "========================================"
echo "Download step complete!"
echo "========================================"
echo

echo "Downloaded genomes:"
ls -lh genome_sequences/*.fa
echo

echo "Original/unfiltered protein files:"
ls -lh unfiltered_protein_sequences/*.fa
echo

echo "Annotations:"
ls -lh annotations/*.gff
echo
