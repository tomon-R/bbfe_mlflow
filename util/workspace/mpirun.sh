gedatsu_simple_mesh_partitioner -n 8
gedatsu_bc_partitioner_R -n 8 -ig node.dat -i D_bc_v.dat
gedatsu_dist_val_partitioner_R -n 8 -ig node.dat -i levelset.dat 
mpirun -np 8 mlflow_fs_sups ./ | tee output.log
