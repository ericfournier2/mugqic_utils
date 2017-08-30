BEDGRAPH=REMOVE
HOMERPREPARSED=REMOVE

# Read command line arguments
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -b|--keepbedgraph)
    BEDGRAPH=KEEP
    ;;
    -h|--keephomerpreparsed)
    HOMERPREPARSED=KEEP
    ;;
esac
shift # past argument or value
done

# Delete intermediaries with no intrinsic value
rm -rf trim

# Delete redundant alignments.
# By default only mdup.hardclip.bam is kept.
for sample in `ls -d alignment/*/ 2>/dev/null` 
do
    # List all readsets and remove them.
    for readset in `ls -d $sample/*/ 2>/dev/null`
    do
        rm -rf $readset
    done
    
    # Remove merged bam files.
    rm -f $sample/*.merged.bam
done

if [ "$BEDGRAPH" = "REMOVE" ]
then 
    rm -f tracks/*/*.bedGraph.gz
    rm -f report/tracks.zip
fi

if [ "$HOMERPREPARSED" = "REMOVE" ]
then 
    rm -rf report/annotation/*/*/preparsed
fi
