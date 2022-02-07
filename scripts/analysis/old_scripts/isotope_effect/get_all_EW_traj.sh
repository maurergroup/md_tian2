for i in `seq 1 10`
do
    cd traj_${i}
    ./../../../scripts/get_EW_traj.py
    cd ..
done
