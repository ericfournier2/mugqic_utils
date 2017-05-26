mugqicdir=$1
genome=$2
rap=$3

# We'll need the UCSC utils to fetch chromosome sizes.
module load mugqic/ucsc/v346

# Fetch the chromosome sizes.
fetchChromSizes $genome > $genome.chrom.sizes

# Loop over all peak files.
for i in $mugqicdir/peak_call/*/*.*Peak
do
    # Determine the sample name.
    narrowRegex='narrowPeak$'
    broadRegex='broadPeak$'
    if [[ $i =~ $narrowRegex ]]
    then
        samplename=`basename $i .narrowPeak`
    elif [[ $i =~ $broadRegex ]]
    then 
        samplename=`basename $i .broadPeak`
    else
        echo 'Not a peak file.'
    fi
    
    mkdir -p $mugqicdir/tracks/PeakBigbed
    mkdir -p $mugqicdir/jobs
    script=$mugqicdir/jobs/$samplename.make_peak_bigbed.sh
    cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A $rap
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=12
bash scripts/generate_bigbed.sh $i $mugqicdir/tracks/PeakBigbed/$samplename.bb $genome.chrom.sizes
EOF
    workdir=`pwd`
    qsub $script -o $script.stdout -e $script.stderr -d $workdir
done
