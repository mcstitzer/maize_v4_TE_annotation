
FILES=($(ls -1 ~/te_reference_fasta/individual_fasta/DTA_*))

### set up a file name for each individual TE fasta
declare -A FILENAMES
FILENAMES=${FILES[*]}
#echo $FILENAME

for i in $FILENAMES
do
### RUN TARGeT
FILE=$(basename "$i")
#echo $FILE
DNATE="${FILE%.*}"
#echo $DNATE
if [ ! -f ${DNATE}.tir.fa.tab ]
then
echo $i
fi
done


### but to get the index to give slurm array
length=${#FILES[@]}
for (( i=0; i<${length}; i++))
do
FILENAME=${FILES[$i]}
FILE=$(basename "$FILENAME")
DNATE="${FILE%.*}"
if [ ! -f ${DNATE}.tir.fa.tab ]
then
echo $i
fi
done
