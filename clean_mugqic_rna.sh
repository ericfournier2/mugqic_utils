# Set argument default values.
DIRECTIONAL=REMOVE
BEDGRAPH=REMOVE

# Read command line arguments
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -d|--keepdirectional)
    DIRECTIONAL=KEEP
    ;;
    -b|--keepbedgraph)
    BEDGRAPH=KEEP
    ;;
esac
shift # past argument or value
done

# Delete intermediaries with no intrinsic value
rm -rf trim
rm -rf reference.Merged
rm -rf alignment_1stPass

# Delete redundant alignments.
# By default only mdup.hardclip.bam is kept.
for sample in `ls -d alignment/*/ 2>/dev/null`
do
    # List all readsets and remove them.
    for readset in `ls -d $sample/*/ 2>/dev/null`
    do
        rm -rf $readset
    done
    
    # Remove everything but the filtered reads
    rm -f $sample/*QueryNameSorted.bam
    rm -f $sample/*.sorted.bam
    rm -f $sample/*.sorted.mdup.bam
    
    if [ $DIRECTIONAL == "REMOVE" ]
    then 
        rm -f $sample/*.forward.bam
        rm -f $sample/*.reverse.bam
    fi
done

# Remove rRNA alignments from the metrics folder.
rm metrics/*/*/*.bam

# Remove bedgraph tracks.
if [ $BEDGRAPH == "REMOVE" ]
then 
    rm tracks/*/*.bedGraph*
    rm tracks.zip
    rm report/tracks.zip
fi
