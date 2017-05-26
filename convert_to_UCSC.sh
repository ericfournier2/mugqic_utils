# Input and output files should be provided on the command line.
INPUTFILE=$1
OUTPUTFILE=$2

# Load required modules.
module load mugqic/samtools
module load mugqic/python/2.7.12

# Script to add chr obtained from https://www.biostars.org/p/10062/
# Add chr prefix to bam file header.
samtools view -H $INPUTFILE |
    perl -lpe 's/Mito/M/' |
    perl -lpe 's/SN:([I]+|I*VI*|I*XV*I*V*|[0-9]+|[XY]|M)\b/SN:chr\$1/' > $INPUTFILE.UCSC.sam

# Add chr prefix to the RNAME and RNEXT columns of alignments.
samtools view $INPUTFILE |
  perl -lpe 's/Mito/M/' |
  perl -lan scripts/fix_chr.pl >> $samplename.sam

# Convert to BAM and index.
samtools view -b -h $samplename.sam > $OUTPUTFILE.bam
samtools index $OUTPUTFILE.bam $OUTPUTFILE.bai

# Remove temporary sam file.
rm $samplename.sam
