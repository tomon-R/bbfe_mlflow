#!/bin/bash

# カラー定義
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
BLUE="\e[34m"
MAGENTA="\e[35m"
RESET="\e[0m"

# 実行と色付け
mlflow_fs_sups | sed "s/^/${RED}/;s/$/${RESET}/" | tee output_1.log &
mlflow_fs_sups | sed "s/^/${GREEN}/;s/$/${RESET}/" | tee output_2.log &
mlflow_fs_sups | sed "s/^/${YELLOW}/;s/$/${RESET}/" | tee output_3.log &
mlflow_fs_sups | sed "s/^/${BLUE}/;s/$/${RESET}/" | tee output_4.log &
mlflow_fs_sups | sed "s/^/${MAGENTA}/;s/$/${RESET}/" | tee output_5.log &

# 全プロセスの終了待ち
wait
