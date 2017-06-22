#!/bin/bash

# Set argument default values.
MUGQIC_DIR=.
INPUTSUFFIX=.sorted.dup.bam
OUTPUTSUFFIX=""
RAP=$RAP_ID
OVERWRITE=FALSE
GENOME=""
HUBNAME="Hub"
HUBLABELSHORT="Hub"
HUBLABELLONG="Hub"
OPTIONFILE="/dev/null"

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
    -w|--overwrite)
    OVERWRITE=TRUE
    ;;
    -n|--hubname)
    HUBNAME=$2
    shift
    ;;
    --hubshortlabel)
    HUBLABELSHORT=$2
    shift
    ;;
    --hublonglabel)
    HUBLABELLONG=$2
    shift
    ;;
    --huboptions)
    OPTIONFILE=$2
    shift
    ;;    
    -g|--genome)
    GENOME=$2
    shift
    ;;
esac
shift # past argument or value
done

# First, generate bigwig files.
bigwigjobs=`bash -x $MUGQIC_UTILS_PATH/generate_bigwig_chip.sh -d $MUGQIC_DIR -s $INPUTSUFFIX -o $OUTPUTSUFFIX -r $RAP -w $OVERWRITE`

# Second, generate bigbed files.
bigbedjobs=`$MUGQIC_UTILS_PATH/generate_peak_bigbeds.sh $MUGQIC_DIR $GENOME $RAP`

# Third, generate the hub proper.
# This will have to be run once all of the above jobs are finished.
mkdir -p $MUGQIC_DIR/jobs
script=$MUGQIC_DIR/jobs/generate_hub.sh
        cat <<EOF > $script
#!/bin/bash
#PBS -N generate_hub.sh
#PBS -A $RAP
#PBS -l walltime=0:15:00
#PBS -l nodes=1:ppn=1

# Load python, which will be required for hub db generation.
module load mugqic/python/2.7.12

# Create hub directory.
mkdir -p $MUGQIC_DIR/hub/$GENOME

# Create symbolic links for bigwigs and bigbeds.
pushd $MUGQIC_DIR/hub/$GENOME
ln -s ../../tracks/*.bw ./
ln -s ../../tracks/PeakBigbed/*.bb ./
popd

# Invoke hub db file generation script.
$MUGQIC_UTILS_PATH/generate_hub_db_files.py -n $HUBNAME -s $HUBLABELSHORT -l $HUBLABELLONG -o $OPTIONFILE $MUGQIC_DIR

EOF

workdir=`pwd`
qsub $MUGQIC_DIR/jobs/generate_hub.sh -o $script.stdout -e $script.stderr -d $workdir -W afterok:$bigwigjobs:$bigbedjobs