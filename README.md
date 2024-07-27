# bbfe_mflow

有限要素法とレベルセット法を使用した二層流れ解析プログラム

## ディレクトリ
二層流ソルバーに関連するディレクトリとファイルは以下の通り。
```plaintext
.
├── bin/
├── docs/
│   └── mlflow_manual         %二層流の計算に関するマニュアル
├── FE_solvers/
│   ├── mlflow_fs             %分離型解法の流体ソルバーを使った二層流解析フォルダ（開発は一体型に途中で移行）
│   └── mlflow_fs_sups        %一体型解法の流体ソルバーを使った二層流解析フォルダ
│   	├── fluid_core.c      %流体の値の更新など
│   	├── fluid_core.h
│   	├── fluid_elemmat.c　 %流体の計算のための行列やベクトルの構築
│   	├── fluid_elemmat.h
│   	├── mlflow_elemmat.c  %レベルセット関数の移流や再初期化のための行列やベクトルの構築
│   	├── mlflow_elemmat.h
│   	├── mlflow_fs_sups.c　%メインの計算部分
│   	├── mlflow_impfunc.c　%レベルセット関数の更新や近似ヘビサイド関数の計算など
│   	├── mlflow_impfunc.h
│   	├── mlflow_utils.c    %ベンチマークと比較のための計測や結果ファイル出力など
│   	└── mlflow_utils.h
├── util/
│   └── workspace
│   	├── data              %ベンチマークの比較データ(csvファイル)
│   	│   ├── 3d-bubble
│   	│   ├── dambreak
│   	│   └── sloshing
│   	├── *.plt             %ベンチマークの結果のグラフ作成
│   	└── *.sh              %ベンチマークの入力ファイル作成シェルスクリプト
├── Makefile
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