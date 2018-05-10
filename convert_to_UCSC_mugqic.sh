#!/bin/bash

# Set argument default values.
MUGQIC_DIR=.
INPUTSUFFIX=.sorted.dup.bam
OUTPUTSUFFIX=".ucsc.bam"
RAP=$RAP_ID
OVERWRITE=FALSE

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
    -o|--overwrite)
    OVERWRITE=TRUE
    ;;
    -l|--slurm)
    SLURM=TRUE
    ;;    
esac
shift # past argument or value
done

for i in $MUGQIC_DIR/alignment/*/*$INPUTSUFFIX
do
    samplename=`basename $i $INPUTSUFFIX`
    sampledir=`dirname $i`
    if [ ! -e $MUGQIC_DIR/tracks/$samplename.bw -o $OVERWRITE -eq TRUE ]
    then     
        mkdir -p $MUGQIC_DIR/jobs
        script=$MUGQIC_DIR/jobs/$samplename.ucsc_convert.sh
        
        cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A $RAP
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
convert_to_UCSC.sh $i $i$OUTPUTSUFFIX
EOF
        workdir=`pwd`
        if [ $SLURM = TRUE ]
        then
            sbatch --time=6:00:00 --account=$RAP -J $script -N 1 --mincpus 16 -o $script.stdout -e $script.stderr -D $workdir
        else
            qsub $script -o $script.stdout -e $script.stderr -d $workdir
        fi
    fi
done
