# Set argument default values.
DIRECTIONAL=REMOVE

# Read command line arguments
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -d|--keepdirectional)
    DIRECTIONAL=KEEP
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
for sample in `ls -d alignment/*/`
do
    # List all readsets and remove them.
    for readset in `ls -d alignment/$sample/*/`
    do
        rm -rf $readset
    done
    
    # Remove everything but the filtered reads
    rm alignment/$sample/*QueryNameSorted.bam
    rm alignment/$sample/*.sorted.bam
    rm alignment/$sample/*.sorted.mdup.bam
    
    if [ $DIRECTIONAL == "REMOVE" ]
    then 
        rm alignment/$sample/*.forward.bam
        rm alignment/$sample/*.reverse.bam
    fi
done
