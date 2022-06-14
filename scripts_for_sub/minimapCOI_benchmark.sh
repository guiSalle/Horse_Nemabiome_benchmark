#!/bin/bash

fwd_primer_C12="RGCHAARCCNGGDYTRTTRYTDGG" #FWD R1_F2
rev_primer_C12="GYTCYAAHGAAATHGAHCTHCTHCG" ## NC2
fwd_primer_rev_C12="CCHARYAAYARHCCNGGYTTDGCY"
rev_primer_rev_C12="CGDAGDAGDTCDATTTCDTTRGARC"

indir=/work2/project/isp/MPN/CYATHOMIX/NMB/mock_raw/
outdir=/work2/project/isp/MPN/CYATHOMIX/NMB/mock_cutadapt/ ### this is produced once (nmb_mock.sh)
wk=/work2/project/isp/MPN/CYATHOMIX/NMB_MOCK/

cd $indir

sampC12=`ls .|grep COI|awk 'BEGIN{FS="_R"}{print $1}'|sort -u|tr "\n" " "`

cd $wk

module load bioinfo/seqkit-v0.16.0
module load bioinfo/seqtk-1.3 

####------ Prep FASTA COI 
# COI
sed 's/ /_/g' ~/save/CYAT/fasta*/Cya_Str_cox1_BOLD_GB_mito.fasta | sed 's/>/>COI/g' > tmp
seqkit rmdup tmp > ./all_cox1_concat.fasta
grep '>' all_cox1_concat.fasta |sed 's/>//g' | awk '{if(index($1,"Quilonia")==0 && index($1,"Murshidia")==0 && index($1,"Khalilia")==0 && index($1,"Kiluluma")==0){print $0}}' > idcox1_noeleph.lst
seqtk subseq all_cox1_concat.fasta idcox1_noeleph.lst > tmp.fa
sed -e 's/AP017681/Cylicostephanus_goldi/g' -e 's/AP017698/Strongylus_vulgaris/g' -e 's/GQ888712/Cylicocyclus_insigne/g' -e 's/GQ888717/Strongylus_vulgaris/g' -e 's/NC_026729/Triodontophorus_brevicauda/g' -e 's/NC_026868/Strongylus_equinus/g' -e 's/NC_031516/Triodontophorus_serratus/g' -e 's/NC_031517/Triodontophorus_nipponicus/g' -e 's/NC_032299/Cylicocyclus_nassatus/g' -e 's/NC_035003/Cyathostomum_catinatum/g' -e 's/NC_035004/Cylicostephanus_minutus/g' -e 's/NC_035005/Poteriostomum_imparidentatum/g' -e 's/NC_038070/Cyathostomum_pateratum/g' -e 's/NC_039643/Cylicocyclus_radiatus/g' -e 's/NC_042141/Cylicodontophorus_bicoronatus/g' -e 's/NC_042234/Coronocyclus_labiatus/g' -e 's/NC_043849/Cylicocyclus_auriculatus/g' -e 's/NC_046711/Cylicocyclus_ashworthi/g' tmp.fa > tmp 
seqkit rmdup tmp > ./all_cox1_concat_clean_rmdup.fasta
rm all_cox1_concat.fasta


####------ Prep FASTA COI 

##rm elephant species
grep '>' all_cox1_concat_clean_rmdup.fasta |sed 's/>//g' | awk '{if(index($1,"Quilonia")==0 && index($1,"Murshidia")==0 && index($1,"Khalilia")==0 && index($1,"Kiluluma")==0){print $0}}' > idcox1_noeleph.lst
seqtk subseq all_cox1_concat_clean_rmdup.fasta idcox1_noeleph.lst > tmp.fa

## Remove duplicate entries (by name)
seqkit rmdup -s < tmp.fa > all_cox1_concat_noeleph_rmdup_clean.fasta

## Remove identical sequences (by sequence)
unclufa=${wk}all_cox1_concat_noeleph_rmdup_clean.fasta

module load bioinfo/usearch11.0.667 ; usearch -fastx_uniques $unclufa -fastaout uniques.fasta

## Cluster db with 100% similarity
for id in 97 99 
do
    usearch -cluster_smallmem uniques.fasta -sortedby other -id 0.$id -centroids all_cox1_concat_noeleph_rmdup_clust$id.fasta
done
id=100
usearch -cluster_smallmem uniques.fasta -sortedby other -id 0.$id -centroids all_cox1_concat_noeleph_rmdup_clust$id.fasta

rm *tmp* all_cox1_concat_noeleph_rmdup_clean.fasta all_cox1_concat_clean_rmdup.fasta


######----- Ready 
cd $wk

echo $sampC12

##------------------- 1. Cutadapt --- done already (nmb_mock.sh)

#for i in $sampC12
#do
#    echo $i
#    sbatch -J"cut.$i" -e e.cut.$i.c12 -o o.cut.$i.c12 --mem=4G -t 00:10:00 --wrap="module load bioinfo/cutadapt-1.14-python-3.4.3; cutadapt --pair-filter any --no-indels \
#-g $fwd_primer_C12 -a $rev_primer_rev_C12 \
#-G $rev_primer_C12 -A $fwd_primer_rev_C12 \
#-m 50 --max-n 1 -q 15 -n 2 --discard-untrimmed \
#-o ${outdir}${i}_R1.fastq.gz -p ${outdir}${i}_R2.fastq.gz ${indir}${i}_R1.fastq.gz ${indir}${i}_R2.fastq.gz"
#done

##---- MINIMAP2 on merged reads
cd $outdir
#gzip *.fastq
cd $wk

for id in 99 97
do
    for B in 1 2 3 4
    do
	for k in 10 13 15
	do
	    for w in 8 9 10
	    do
		fastaCOI=${wk}all_cox1_concat_noeleph_rmdup_clust$id.fasta
		wk2=${wk}mock_MINIMAPm_clus${id}_B${B}_k${k}_w${w}
		mkdir $wk2
		cd $wk2
		echo $id $B $k $w $wk2
		### No filtering for mapping quality but keep proper pairs only/ filtering MQ15
		for sample in $sampC12
		do
		    reffa=$fastaCOI
#		    sbatch -J"minimap" -o mini.$sample.o -e mini.$sample.e -c 4 --mem=10G -t 00:15:00 --wrap="module load bioinfo/usearch11.0.667; module load bioinfo/minimap2-2.11; module load bioinfo/samtools-1.10; cd $outdir; usearch -fastq_mergepairs ./${sample}_R1.fastq -relabel @ -fastqout ./$sample.merged.fq ; cd ${wk2}; minimap2 -t 4 -ax sr -A2 -B$B -O4,24 -E2,1 -k$k -w$w $reffa ${outdir}${sample}.merged.fq | samtools view -b -F 2308 -q 0| samtools sort -@ 4 - -o ${sample}_sorted.bam; samtools flagstat ${sample}_sorted.bam > ${sample}_sorted.bam.flagstat.txt ; samtools index ${sample}_sorted.bam ; samtools view -b -F 2308 -q 30 ${sample}_sorted.bam -o ${sample}_filtered_sorted.bam ; samtools index ${sample}_filtered_sorted.bam ; module load bioinfo/bedtools2-2.29.0 ; bedtools genomecov -d -ibam ${sample}_filtered_sorted.bam | gzip > ${sample}_genomecov_d.txt.gz ; bedtools genomecov -d -ibam ${sample}_sorted.bam | gzip > ${sample}_nofilt_genomecov_d.txt.gz"
		done
	    done
	done
    done
done

### Compute # of reads per sequence for the chosen set of parameters and full range evaluation on community composition prediction
id=99
B=4
k=10
w=10
for sample in $sampC12
do
    cd ${wk}mock_MINIMAPm_clus${id}_B${B}_k${k}_w${w}
    echo $sample
    module load bioinfo/samtools-1.11 ; samtools idxstats ${sample}_sorted.bam > ${sample}_idxstats.txt
done
