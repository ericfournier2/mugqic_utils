inputbed=$1
outputbb=$2
chromsizes=$3

module load mugqic/ucsc

cut -f 1-4 $inputbed | sort -k1,1 -k2,2n  > $inputbed.sorted
bedToBigBed $inputbed.sorted $chromsizes $outputbb
rm $inputbed.sorted