for i in `seq 1 10`
do
	cd traj_$i
	nohup ../../../md_tian/md_tian2.x md_tian.inp > output.dat &
	cd ..
done
