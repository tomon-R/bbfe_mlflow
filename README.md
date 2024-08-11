# bbfe_mflow

有限要素法とレベルセット法を使用した二層流れ解析プログラム

## ディレクトリ
二層流ソルバーに関連するディレクトリとファイルは以下の通り。
```plaintext
.
├── bin/
├── docs/
│   └── mlflow_manual          # 二層流の計算に関するマニュアル
├── FE_solvers/                # 解析ソルバーのディレクトリ
│   ├── mlflow_fs              # 分離型解法の流体ソルバーを使った二層流解析ディレクトリ（開発は一体型に途中で移行）
│   └── mlflow_fs_sups         # 直接型解法の流体ソルバーを使った二層流解析ディレクトリ
│   	├── fluid_core.c       # 流体の値の更新など
│   	├── fluid_core.h
│   	├── fluid_elemmat.c    # 流体の計算のための行列やベクトルの構築
│   	├── fluid_elemmat.h
│   	├── mlflow_elemmat.c   # レベルセット関数の移流や再初期化のための行列やベクトルの構築
│   	├── mlflow_elemmat.h
│   	├── mlflow_fs_sups.c   # 二層流解析のmain関数などのメインの計算部分
│   	├── mlflow_impfunc.c   # レベルセット関数の更新や近似ヘビサイド関数の計算など
│   	├── mlflow_impfunc.h
│   	├── mlflow_utils.c     # ベンチマークと比較のための計測や結果ファイル出力など
│   	├── mlflow_utils.h
│   	└── Makefile
├── util/
│   └── workspace              # 入力ファイル作成などの作業用ディレクトリ
│   	├── mesh
│   	│   └── levelset_gen.c # レベルセット関数の初期状態の入力ファイルlevelset.datの作成プログラム
│   	├── data               # ベンチマークの比較データ(csvファイル)
│   	│   ├── 3d-bubble
│   	│   ├── dambreak
│   	│   └── sloshing
│   	├── *.plt              # ベンチマークの結果のグラフ作成
│   	└── *.sh               # ベンチマークの入力ファイル作成シェルスクリプト
└── README.md
```

## インストール・コンパイル方法
### 前提
以下のパッケージが入っていること。

- make
- cmake
- git
- gcc (gfortran)
- MPI

### インストール手順
```shell
$ monolis_install.sh
$ install.sh
$ util_install.sh
```

### コンパイル
```shell
$ cd ./FE_solvers/mlflow_fs_sups
$ make
```

## 実行方法
### 入力ファイルの準備
```shell
$ cd ./util/workspace #入力ファイル作成などの作業用ディレクトリに移動
$ make_sloshing_allnoslip.sh 40 8 48 1.0 0.2 1.2 #一例としてスロッシングの入力ファイルを造るシェルの実行
$ ls #Sloshing_40_8_1.0_0.2_1.2のようなフォルダが作成され、中にnode.dat, elem.dat, levelset.dat, D_bc_v.datがあることを確認する
$ cp Sloshing_40_8_1.0_0.2_1.2/*.dat ../FE_solvers/mlflow_fs_sups/sloshing #mlflow_fs_sups下に解析用フォルダ(sloshing)を作っておきそこに入れる
$ #cond.datファイルは最初から準備されているか、なければ自分で作成し解析用フォルダに入れておく
```

### 解析条件入力ファイル(cond.dat)

```plaintext
#num_ip_each_axis 1   # 数値積分の1辺の積分点数
2                        
#mat_epsilon 1        # 線形ソルバーの反復法の許容誤差
1.000000000000000e-08    
#mat_max_iter 1       # 線形ソルバーの最大反復回数
10000                    
#time_spacing 1       # 流体ソルバーの時間刻み幅
0.001                    
#finish_time 1        # 計算終了時間
10                       
#output_interval 1    # 結果ファイルを出力するステップ間隔
100
#density_l 1          # 液相の密度
1000
#density_g 1          # 気相の密度
1.2
#viscosity_l 1        # 液相の粘性係数
1.0e-3
#viscosity_g 1        # 気相の粘性係数
1.8e-5
#gravity 3            # 重力加速度ベクトル
0.0
0.0
-9.8
#size_interface 1     # 界面厚さ
0.1
#surf_tension_coef 1  # 表面張力係数
0
#dt_reinit 1          # レベルセット関数の再初期化計算の時間刻み幅dt
1e-4
#epsilon_reinit 1     # レベルセット関数の再初期化計算のパラメータ（sgn(phi)を近似するパラメータ）
0.1
#delta_reinit 1       # レベルセット関数の再初期化計算の収束判定値
1e-4
#max_iter_reinit      # レベルセット関数の再初期化の最大反復回数
5
#accel_amp 3          # 慣性加速度ベクトルの振幅
0.0093
0
0
#accel_angle_vel 3    # 慣性加速度ベクトルの角速度
5.311
0
0
#output_option 1      # ベンチマークの計測ファイルを出力するオプション（1:ダムブレイク、2:気泡上昇流れ、3:スロッシング）
3
#ale_option 1         # ALE法でメッシュを移動させるかどうかのオプション (0:メッシュ移動させない、1:メッシュ移動させる)
0
```

### 解析の実行
```shell
$ cd ../FE_solvers/mlflow_fs_sups/ #ソルバーのディレクトリに移動
$ ./mlflow_fs_sups ./sloshing #解析用フォルダ(sloshing)に結果ファイルが保存されていく
```

## 解析結果の確認
### 可視化

- 解析用フォルダに.vtkファイルが出力されるのでparaviewなどで可視化できる

### グラフ化

- 事前にgnuplotを公式ホームページからインストールする(http://www.gnuplot.info/)
- 解析用フォルダにcond.datの#output_typeの番号によってはベンチマークとの比較のための計測量のcsvが出力される
- util/workspace/data/sloshingなどのフォルダにcsvファイルを入れる
- util/workspace/graph_sloshing.pltを実行することでグラフ化できる

## ベンチマーク
### ダムブレイク

### 気泡上昇流れ

### スロッシング