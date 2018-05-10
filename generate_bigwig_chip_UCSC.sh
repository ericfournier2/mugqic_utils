for i in output/pipeline/alignment/*/*.sorted.dup.bam
do
    samplename=`basename $i .sorted.dup.bam`
    if [ ! -e output/pipeline/tracks/$samplename.bw ]
    then
       mkdir -p output/pipeline/jobs
       script=output/pipeline/jobs/$samplename.make_bigwig_ucsc.sh
       cat <<EOF > $script
#!/bin/bash
#PBS -N $script
#PBS -A eav-760-aa
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=16
module load mugqic/samtools
module load mugqic/python/2.7.12
# Script to add chr obtained from https://www.biostars.org/p/10062/
# Add chr prefix to bam file.
samtools view -H $i |
    perl -lpe 's/Mito/M/' |
    perl -lpe 's/SN:([I]+|I*VI*|I*XV*I*V*|[0-9]+|[XY]|M)\b/SN:chr\$1/' > $samplename.sam

samtools view $i |
  perl -lpe 's/Mito/M/' |
  perl -lan scripts/fix_chr.pl >> $samplename.sam
samtools view -b -h $samplename.sam > $samplename.bam
samtools index $samplename.bam $samplename.bai
rm $samplename.sam
bamCoverage -e 200 --binSize 5 -p 16 --normalizeUsingRPKM \
    -b $samplename.bam \
    -o output/pipeline/tracks/$samplename.bw
rm $samplename.bam
rm $samplename.bai
EOF
        workdir=`pwd`
        qsub $script -o $script.stdout -e $script.stderr -d $workdir
    fi
done
