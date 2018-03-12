# Set default values.
MUGQIC_DIR="."
GENOME="hg19"
RAP=$RAP_ID
MUGQIC_UTILS_HOME=~/mugqic_utils
DRY=FALSE

# Read command line arguments
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -d|--mugqicdir)
    MUGQIC_DIR=$2
    shift
    ;;
    -s|--inputsuffix)
    INPUTSUFFIX=$2
    shift
    ;;
    -o|--outputsuffix)
    OUTPUTSUFFIX=$2
    shift
    ;;    
    -r|--rapid)
    RAP=$2
    shift
    ;;
    -g|--genome)
    GENOME=$2
    shift
    ;;    
    -w|--overwrite)
    OVERWRITE=TRUE
    ;;
    -l|--slurm)
    SLURM=TRUE
    ;;
    -y|--dry)
    DRY=TRUE
    ;;    
esac
shift # past argument or value
done

# We'll need the UCSC utils to fetch chromosome sizes.
module load mugqic/ucsc

# Fetch the chromosome sizes.
fetchChromSizes $GENOME > $GENOME.chrom.sizes

# Loop over all peak files.
jobids=""
for i in $MUGQIC_DIR/peak_call/*/*.*Peak
do
    # Determine the sample name.
    narrowRegex='narrowPeak$'
    broadRegex='broadPeak$'
    samplename="ERROR"
    if [[ $i =~ $narrowRegex ]]
    then
        samplename=`basename $i .narrowPeak`
    elif [[ $i =~ $broadRegex ]]
    then 
        samplename=`basename $i .broadPeak`
    else
        echo 'Not a peak file.'
    fi
    
    mkdir -p $MUGQIC_DIR/tracks/PeakBigbed
    mkdir -p $MUGQIC_DIR/jobs

    if [ "$samplename" != "ERROR" ]
    then 
        script=$MUGQIC_DIR/jobs/$samplename.make_peak_bigbed.sh
        cat <<EOF > $script
#!/bin/bash

bash $MUGQIC_UTILS_HOME/generate_bigbed.sh $i $MUGQIC_DIR/tracks/PeakBigbed/$samplename.bb $GENOME.chrom.sizes
EOF
        workdir=`pwd`
        if [ $DRY = FALSE ]
        then
            if [ $SLURM = TRUE ]
            then
                jobid=`sbatch --time=6:00:00 --account=$RAP -J $script -N 1 --mem=32G --mincpus 16 -o $script.stdout -e $script.stderr -D $workdir $script `
                jobid=`echo $jobid | sed -e 's/Submitted batch job //'`
                jobids=$jobids:$jobid
            else
                jobids=$jobids:`qsub $script -o $script.stdout -e $script.stderr -d $workdir -N $script -A $RAP -l walltime=6:00:00 -l nodes=1:ppn=12`
            fi    
        fi
    fi
done

echo $jobids
