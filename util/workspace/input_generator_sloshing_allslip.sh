#!/bin/bash

if [ $# -ne 6 ];then
    echo "Usage: num_elem_x num_elem_y num_elem_z len_x len_y len_z";
    exit 1
fi

../meshgen/meshgen_hex $1 $2 $3 $4 $5 $6

../surface/surf_conn 1

../mesh/mesh_surf_extract 0.0 0.0 0.0 0.0 $5  $6  -oe surf_left.dat -ov surf_left.vtk
../mesh/mesh_surf_extract $4  0.0 0.0 $4  $5  $6  -oe surf_right.dat -ov surf_right.vtk
../mesh/mesh_surf_extract 0.0 0.0 0.0 $4  0.0 $6  -oe surf_front.dat -ov surf_front.vtk
../mesh/mesh_surf_extract 0.0 $5  0.0 $4  $5  $6  -oe surf_back.dat -ov surf_back.vtk
../mesh/mesh_surf_extract 0.0 0.0 0.0 $4  $5  0.0 -oe surf_bottom.dat -ov surf_bottom.vtk
../mesh/mesh_surf_extract 0.0 0.0 $6  $4  $5  $6  -oe surf_top.dat -ov surf_top.vtk

../surface/surf_dbc_sups 4 0.0 -1 -1 -1 -ie surf_left.dat  -o D_bc_vx1.dat
../surface/surf_dbc_sups 4 0.0 -1 -1 -1 -ie surf_right.dat  -o D_bc_vx2.dat
../surface/surf_dbc_sups 4 -1 0.0 -1 -1 -ie surf_front.dat  -o D_bc_vy1.dat
../surface/surf_dbc_sups 4 -1 0.0 -1 -1 -ie surf_back.dat  -o D_bc_vy2.dat
../surface/surf_dbc_sups 4 -1 -1 0.0 -1 -ie surf_bottom.dat  -o D_bc_vz1.dat
../surface/surf_dbc_sups 4 -1 -1 0.0 -1 -ie surf_top.dat  -o D_bc_vz2.dat

../surface/surf_bc_merge   D_bc_vx1.dat D_bc_vx2.dat  -o D_bc_vx.dat
../surface/surf_bc_merge   D_bc_vy1.dat D_bc_vy2.dat  -o D_bc_vy.dat
../surface/surf_bc_merge   D_bc_vz1.dat D_bc_vz2.dat  -o D_bc_vz.dat
../surface/surf_bc_merge   D_bc_vx.dat D_bc_vz.dat  -o D_bc_vxz.dat
../surface/surf_bc_merge   D_bc_vxz.dat D_bc_vy.dat  -o D_bc_v.dat

OPTION=1
x_min=0
x_max=$(echo "scale=5; $4 * 2" | bc)
y_min=0
y_max=$5
z_min=0
z_max=$(echo "scale=5; 0.5" | bc)
#z_max=$(echo "scale=5; $6 / 2.0" | bc)

echo -e "../mesh/levelset_gen $OPTION $x_min $y_min $z_min $x_max $y_max $z_max \n"
../mesh/levelset_gen $OPTION $x_min $y_min $z_min $x_max $y_max $z_max

DIR_NAME="Sloshing_"$1"_"$2"_"$3"_"$4"_"$5"_"$6

echo -e " mkdir $DIR_NAME \n"
mkdir $DIR_NAME

mv ./*.dat $DIR_NAME
mv ./*.vtk $DIR_NAME

rm $DIR_NAME/D_bc_vx*.dat
rm $DIR_NAME/D_bc_vy*.dat
rm $DIR_NAME/D_bc_vz*.dat
rm $DIR_NAME/surf_*.dat