for i in 1 2 3 4 5 6 7
do
for j in 1 2 3 4
do
for k in 1 2 3
do
for z in 1 2 3
do
qsub -cwd run.sh $i $j $k $z
done
done
done
done