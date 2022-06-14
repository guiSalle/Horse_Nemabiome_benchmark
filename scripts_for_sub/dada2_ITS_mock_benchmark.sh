###------------------------------------------------
###- gs12 --- CYATHOMIX barcode data ITS-2 and COI
###------------------------------------------------

bc='ITS'

##------------------- 1. Cutadapt
#fwd_primer_ITS="ACGTCTGGTTCAGGGTTGTT" ## NC1
#rev_primer_ITS="TTAGTTTCTTTTCCTCCGCT" ## NC2
#fwd_primer_rev_ITS="AACAACCCTGAACCAGACGT"
#rev_primer_rev_ITS="AGCGGAGGAAAAGAAACTAA"

#indir=/work2/project/isp/MPN/CYATHOMIX/NMB_XP/raw/
#outdir=/work2/project/isp/MPN/CYATHOMIX/NMB_XP/cutadapt/
#wk=/work2/project/isp/MPN/CYATHOMIX/NMB_XP/

#cd $indir

#sampITS=`ls .|grep ITS|awk 'BEGIN{FS="_R"}{print $1}'|sort -u|tr "\n" " "`

#cd $wk
#for i in $sampITS
#do
#    echo $i
#    sbatch -J"cut.$i" -e e.cut.$i.its -o o.cut.$i.its --mem=4G -t 00:10:00 --wrap="module load bioinfo/cutadapt-1.14-python-3.4.3; cutadapt --pair-filter any --no-indels \
#-g $fwd_primer_ITS -a $rev_primer_rev_ITS \
#-G $rev_primer_ITS -A $fwd_primer_rev_ITS \
#-m 50 --max-n 1 -q 15 -n 2 --discard-untrimmed \
#-o ${outdir}${i}_R1.fastq.gz -p ${outdir}${i}_R2.fastq.gz ${indir}${i}_R1.fastq.gz ${indir}${i}_R2.fastq.gz"
#done

####---- Directories for subsequent analyses
indir=/work/project/isp/MPN/CYATHOMIX/NMB_XP/cutadapt/ ### Prepared with nmb_cya_sainf.sh
work=/work/project/isp/MPN/CYATHOMIX/NMB_XP/

maniftot=/work/project/isp/MPN/CYATHOMIX/NMB_XP/all_DADA2/manifest.csv
manif=${work}manifest$bc.csv

###---------- ITS-2 data for life-stages
#cd $indir
#ls *$bc*.fastq.gz|awk -v d=$indir 'BEGIN{FS="_"}{print $1","d""$0",direction"substr($4,1,2)}'|sed 's/directionR1/forward/g'|sed 's/directionR2/reverse/g' >> $maniftot

echo "sample-id,absolute-filepath,direction" > $manif
grep ITS $maniftot |awk 'BEGIN{FS=","}{if(index($1,"ITSM")!=0 || length(substr($1,index($1,"-"))+1) > 3 || $1=="ITSH2O"){print $0}}' >> $manif

cd $work

### Poissant et al. used BAND of -1 (full WNeedleman-Wunsch alignment)
### For ITS-2 32 is recommended // 16 is the default
### We use 32 as it performs best for benchmarking

input_Fdir=${work}indirF/
input_Rdir=${work}indirR/

## Create sym links for files
mkdir $input_Fdir
cd $input_Fdir
for i in `awk 'BEGIN{FS=","}{if(NR>1){print $2}}' $manif |grep R1.|tr "\n" " "`
do
    ln -s $i .
done
mkdir $input_Rdir
cd $input_Rdir
for i in `awk 'BEGIN{FS=","}{if(NR>1){print $2}}' $manif |grep R2.|tr "\n" " "`
do
    ln -s $i .
done

cd $work
nt=20

for BAND_SIZE in 32 #-1 16 32
do
    for truncF in 200 #217
    do
	mxeeF=1
	mxeeR=1
	outdir=${work}dada2_R_mxee${mxeeF}${mxeeR}_trunc${truncF}_BS$BAND_SIZE/
	mkdir $outdir
	filtered_Fdir=${outdir}data_filtF
	mkdir $filtered_Fdir
	filtered_Rdir=${outdir}data_filtR
	mkdir $filtered_Rdir
	outtsv=${outdir}output.tsv
	tracktsv=${outdir}track.tsv
	trimLeft=0

	echo "$input_dir $outtsv $tracktsv $filtered_dir $nt $BAND_SIZE"
#	sbatch -c $nt --mem=20G -e e.dada2.$truncF.$mxeeF.$BAND_SIZE -o o.dada2.$truncF.$mxeeR.$BAND_SIZE -J"BS$BAND_SIZE" --wrap="module load compiler/gcc-7.2.0 ; module load system/R-4.0.2_gcc-7.2.0 ; Rscript ~/save/scripts/run_dada_paired.R $input_Fdir $input_Rdir $outtsv $tracktsv $filtered_Fdir $filtered_Rdir $truncF $truncF 0 0 $mxeeF $mxeeR 2 12 pseudo consensus 1.0 $nt 1000000 $BAND_SIZE ; sed 's/#OTU ID/OTUID/g' $outtsv > tmp ; mv tmp $outtsv"
    done
done


