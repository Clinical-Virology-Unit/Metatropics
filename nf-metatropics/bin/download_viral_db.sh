#!/bin/bash

# Create database directory
mkdir -p ViralRefseq
cd ViralRefseq

# Download viral sequences from NCBI
echo "Downloading viral sequences from NCBI..."
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz

# Decompress file
gunzip viral.1.1.genomic.fna.gz

# Filter for eukaryotic viruses (excluding phages) and create ref2taxid
echo "Filtering for eukaryotic viruses and creating ref2taxid..."
awk '
BEGIN { RS = ">" }
!/phage/ && !/bacteriophage/ && /Eukaryota/ {
    if (NR > 1) {  # Skip empty first record
        # Extract accession and taxid from header
        acc = $1
        gsub(/\.[0-9]+$/, "", acc)  # Remove version number
        taxid = $0
        sub(/.*taxid=/, "", taxid)
        sub(/ .*/, "", taxid)
        print acc "\t" taxid
        # Print sequence with header
        print ">" $0
    }
}
' viral.1.1.genomic.fna > DB.fa

# Download taxonomy database
echo "Downloading taxonomy database..."
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

# Clean up intermediate files
rm viral.1.1.genomic.fna taxdump.tar.gz

echo "Database preparation complete!"
echo "Files created:"
echo "- DB.fa: Reference sequences"
echo "- Taxonomy files: nodes.dmp, names.dmp, etc." 