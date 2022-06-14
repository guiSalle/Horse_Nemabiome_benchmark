###----------------------------
###- gs12 --- dada2 pipeline
###----------------------------

bc=$1
mxeeF=$2
mxeeR=$3

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters"
fi

indir=/work/project/isp/MPN/CYATHOMIX/NMB/mock_cutadapt/
work=/work/project/isp/MPN/CYATHOMIX/NMB_MOCK/

cd $work

### Poissant et al. used BAND of -1 (full WNeedleman-Wunsch alignment)
### For ITS-2 32 is recommended
### BS is set to 16 by default

input_Fdir=${work}mock_${bc}_indirF/
input_Rdir=${work}mock_${bc}_indirR/

nt=20

for BAND_SIZE in -1 16 32
do
    for truncF in 200 217
    do

	outdir=${work}mock_dada2_R_${bc}_mxee${mxeeF}${mxeeR}_trunc${truncF}_BS$BAND_SIZE/
	mkdir $outdir
	filtered_Fdir=${outdir}data_filtF
	mkdir $filtered_Fdir
	filtered_Rdir=${outdir}data_filtR
	mkdir $filtered_Rdir
	outtsv=${outdir}output.tsv
	tracktsv=${outdir}track.tsv
	trimLeft=0

	echo "$input_dir $outtsv $tracktsv $filtered_dir $nt $BAND_SIZE"
	sbatch -c $nt --mem=20G -e e.dada2.$bc.$truncF.$mxeeF.$BAND_SIZE.%J -o o.dada2.$bc.$truncF.$mxeeR.$BAND_SIZE.%J -J"$bc.$truncF.$mxeeR.$BAND_SIZE" --wrap="module load compiler/gcc-7.2.0 ; module load system/R-4.0.2_gcc-7.2.0 ; Rscript ~/save/scripts/run_dada_paired.R $input_Fdir $input_Rdir $outtsv $tracktsv $filtered_Fdir $filtered_Rdir $truncF $truncF 0 0 $mxeeF $mxeeR 2 12 pseudo consensus 1.0 $nt 1000000 $BAND_SIZE ; sed 's/#OTU ID/OTUID/g' $outtsv > tmp ; mv tmp $outtsv"
    done
done


