# read in file of samples IDs
# for each sample
# 	print father
# 	print mother
# 	print 0
# 	print 0

sample_IDs=$1
family_ID=$2

readarray -r all_samples < $sample_IDs

for sample in "${all_samples[@]}"
do
	parents=( $(biscuit -D PN=$sample biscuit_parent.sh) )
	father=${parents[0]}
	mother=${parents[1]}

	echo $family_ID $sample $father $mother 0 0
done
