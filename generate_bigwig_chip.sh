#!/bin/bash

# Set argument default values.
MUGQIC_DIR=.
INPUTSUFFIX=.sorted.dup.bam
OUTPUTSUFFIX=""
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
esac
shift # past argument or value
done

for i in $MUGQIC_DIR/alignment/*/*$INPUTSUFFIX
do
    samplename=`basename $i $INPUTSUFFIX`
    if [ ! -e $MUGQIC_DIR/tracks/$samplename.bw -o "$OVERWRITE" = "TRUE" ]
    then     
        mkdir -p $MUGQIC_DIR/jobs
        script=$MUGQIC_DIR/jobs/$samplename.make_bigwig.sh
        
        cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A $RAP
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
module load mugqic/python/2.7.12
bamCoverage -e 200 --binSize 5 -p 16 --normalizeUsingRPKM \
    -b $i \
    -o $MUGQIC_DIR/tracks/$samplename$OUTPUTSUFFIX.bw
EOF
        workdir=`pwd`
        qsub $script -o $script.stdout -e $script.stderr -d $workdir
    fi
done
