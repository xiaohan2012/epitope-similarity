#The batch processing file that calculates the pairwise similarity scores for a given set of structures (PDB files)
#
#the following three variables need to be specified by the user:
#
# PDB_DIR: where the pdb files are
# FP_DIR: where to store the finger print directories
# SCORE_DIR: where to store the scores

PDB_DIR=test/data
FP_DIR=fp
SCORE_DIR=score

# HORIZONTAL=(2 4 5)
# VERTICAL=(2 4 5)
# CUTOFF=(20 5 10 15)

HORIZONTAL=(2 4) 
VERTICAL=(2 4)   
CUTOFF=(20 15)

if [ ! -d $SCORE_DIR ]; then
    mkdir -p $SCORE_DIR
else
    rm $SCORE_DIR/* #remove anything new
fi



files=$(ls $PDB_DIR/*)

is_computed () {
    r1=$(cat $1 | grep "$2 $3" | wc -l) 
    r2=$(cat $1 | grep "$3 $2" | wc -l) 

    if [ $r1 == 1 ] || [ $r2 == 1 ]; then
	return 0
    else
	return 1
    fi
}

#generate the fingerprints
echo "Generating fingerprints..."

for h in ${HORIZONTAL[@]}; do
    for v in ${VERTICAL[@]}; do
	fp_dir=$FP_DIR/h$h-v$v
	
	if [ ! -d "$fp_dir" ]; then #if dir not exist, then make it
	    mkdir -p $fp_dir
	fi
	
	for file in $files; do
	    pdb_id=$(basename $file)
	    #python get_fp.py $file $h $v > $fp_dir/$pdb_id.fp
	done
    done
done

echo "Calculating pairwise scores..."
for h in ${HORIZONTAL[@]}; do
    for v in ${VERTICAL[@]}; do
	
	fp_dir=$FP_DIR/h$h-v$v
	
	for cutoff in ${CUTOFF[@]}; do
	    pairwise_score_file="$SCORE_DIR/h$h-v$v-c$cutoff.score"
	    touch $pairwise_score_file
	    
	    for query_pdb in $files; do	    
		query_pdb_id=$(basename $query_pdb)
		query_fp=$fp_dir/$query_pdb_id.fp
		
		for against_pdb in $files; do
		    against_pdb_id=$(basename $against_pdb)
		    against_fp=$fp_dir/$against_pdb_id.fp

		    echo "$query_pdb_id VS $against_pdb_id($h $v $cutoff)"
		    is_computed $pairwise_score_file $query_pdb_id $against_pdb_id
		    if [ $?  == 0 ]; then #we already done that
			echo "skip"
			continue
		    fi
			
		    result=$(python get_sim_score.py --query-pdb=$query_pdb  --against-pdb=$against_pdb  --query-fp=$query_fp --against-fp=$against_fp --cutoff=$cutoff)
		    echo "$query_pdb_id $against_pdb_id $result" >> $pairwise_score_file
		done
	    done
	done
    done
done				
