#!/bin/bash

# ファイルまたはディレクトリのパスを入力
FILE1="../example"
FILE2="../Dambreak_allslip_40_6_40_0.584_0.0876_0.584"

filenames=("cond.dat" "D_bc_v.dat" "elem.dat" "node.dat" "levelset.dat")

for filename in ${filenames[@]}; do
    diff -s "$FILE1/$filename" "$FILE2/$filename" > "$filename.txt"
done