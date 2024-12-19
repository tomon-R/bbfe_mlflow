# BBFEフレームワーク概要

BBFE（Bebops Based Finite Element）は、有限要素法をC言語で実装するためのフレームワークです。

## 目次

- [BBFEフレームワーク概要](#bbfeフレームワーク概要)
  - [目次](#目次)
  - [構造体](#構造体)
  - [入力ファイル](#入力ファイル)
  - [関数ライブラリ](#関数ライブラリ)
    - [**FE\_elemmat**](#fe_elemmat)
    - [**FE\_manusol**](#fe_manusol)
    - [**FE\_std**](#fe_std)
    - [**FE\_sys**](#fe_sys)
    - [**libBB**](#libbb)
  - [メッシュと境界条件の生成](#メッシュと境界条件の生成)
      - [1. メッシュ生成](#1-メッシュ生成)
      - [2. 表面メッシュ抽出](#2-表面メッシュ抽出)
      - [3. 表面の切り出し・除外](#3-表面の切り出し除外)
      - [4. 境界条件の設定](#4-境界条件の設定)
      - [5. 使用例](#5-使用例)
  - [Makefileによるコンパイル](#makefileによるコンパイル)
      - [1. Makefileの概要](#1-makefileの概要)
      - [2. Makefileの変更点](#2-makefileの変更点)
      - [3. コンパイルの手順](#3-コンパイルの手順)
      - [4. `solve_mat`の使用とMakefile](#4-solve_matの使用とmakefile)

***

## 構造体

BBFEフレームワークでは、様々な計算を管理するために複数の構造体が定義されています。主な構造体は以下の通りです。

* **FE\_SYSTEM**: フレームワークの基幹となる構造体。解析に必要な情報を統合し管理します。以下の構造体を含みます。
  * **BBFE\_BASIS**: 形状関数や積分点に関する情報。
  * **BBFE\_DATA**: 節点座標やコネクティビティ、ヤコビアンなど物理空間の情報。
  * **BBFE\_BC**: 境界条件に関する情報。
  * **VALUES**: 材料特性や物理量の値を管理する構造体。解析によってこの構造体を変更することが主な操作となります。
  * **CONDITIONS**: 計算条件（例：積分点の数や境界条件のファイル名など）を管理。
  * **MONOLIS**: 行列計算に用いられる並列計算対応の構造体。

## 入力ファイル

BBFEでは、解析の設定に必要な各種情報を入力ファイルから読み込みます。主な入力ファイルは以下の通りです。

* **node.dat**: 節点データ（メッシュの節点座標を定義）。
* **elem.dat**: 要素データ（要素の節点接続情報を定義）。
* **D\_bc\_v.dat**: 速度のDirichlet境界条件を定義。
* **D\_bc\_p.dat**: 圧力のDirichlet境界条件を定義。
* **cond.dat**: 計算条件を定義。例として、#num\_ip\_each\_axis（積分点数）や#density（密度）などを設定可能です。

## 関数ライブラリ

BBFEは以下の関数ライブラリを使用して、各種有限要素計算を効率的に行います。

***

### **FE\_elemmat**

* **convdiff.c / convdiff.h**
  * `BBFE_elemmat_convdiff_mat_conv` - 移流項の要素剛性行列を計算。
  * `BBFE_elemmat_convdiff_mat_diff` - 拡散項の要素剛性行列を計算。
  * `BBFE_elemmat_convdiff_vec_source` - ソース項（右辺ベクトル）を計算。
  * `BBFE_elemmat_convdiff_stab_coef` - 安定化項の係数を計算。
  * `BBFE_elemmat_convdiff_mat_stab_conv` - 安定化移流項の剛性行列を計算。
  * `BBFE_elemmat_convdiff_vec_stab_source` - 安定化ソース項を計算。
  * `BBFE_elemmat_convdiff_mat_mass` - 質量行列を計算。
  * `BBFE_elemmat_convdiff_vec_mass` - 質量項の右辺ベクトルを計算。
  * `BBFE_elemmat_convdiff_stab_coef_ns` - Navier-Stokes用の安定化係数を計算。
  * `BBFE_elemmat_convdiff_mat_stab_mass` - 質量行列の安定化項を計算。
  * `BBFE_elemmat_convdiff_vec_stab_mass` - 質量項の安定化右辺ベクトルを計算。

* **equivval.c / equivval.h**
  * `BBFE_elemmat_equivval_volume_smooth_function` - 平滑化関数を用いた体積積分を計算。
  * `BBFE_elemmat_equivval_relative_L2_error_scalar` - 相対L2誤差を計算。

* **fluid.c / fluid.h**
  * `BBFE_elemmat_fluid_supg_coef` - SUPG（流線方向アップウィンド）法の係数を計算。
  * `BBFE_elemmat_fluid_sups_coef` - SUPS（流線方向と圧力方向安定化）法の係数を計算。
  * `BBFE_elemmat_fluid_fs_mat_pred_expl` - Fractional step法の予測剛性行列を計算。
  * `BBFE_elemmat_fluid_fs_vec_pred_expl` - Fractional step法の予測右辺ベクトルを計算。
  * `BBFE_elemmat_fluid_fs_vec_ppe` - 圧力ポアソン方程式用右辺ベクトルを計算。
  * `BBFE_elemmat_fluid_fs_vec_corr` - 速度補正の右辺ベクトルを計算。
  * `BBFE_elemmat_fluid_sups_mat` - SUPS法の剛性行列を計算。
  * `BBFE_elemmat_fluid_sups_vec` - SUPS法の右辺ベクトルを計算。

* **set.c / set.h**
  * `BBFE_elemmat_set_Jacobian_array` - 要素のヤコビアン行列を設定。
  * `BBFE_elemmat_set_local_array_scalar` - スカラー場の局所配列を設定。
  * `BBFE_elemmat_set_local_array_vector` - ベクトル場の局所配列を設定。
  * `BBFE_elemmat_set_Jacobi_mat` - ヤコビ行列を設定。
  * `BBFE_elemmat_set_shapefunc_derivative` - 形状関数の微分を設定。
  * `BBFE_elemmat_set_global_mat_cmass_const` - 収束した質量行列を設定。
  * `BBFE_elemmat_set_global_mat_cmass_const_C` - コンスタント質量行列を設定（C用）。
  * `BBFE_elemmat_set_global_mat_Laplacian_const` - コンスタントラプラシアン行列を設定。
  * `BBFE_elemmat_set_global_mat_Laplacian_const_C` - コンスタントラプラシアン行列を設定（C用）。

* **solid.c / solid.h**
  * `BBFE_elemmat_solid_mat_dispstr_linear` - 線形弾性の変位-ひずみ行列を計算。
  * `BBFE_elemmat_solid_mat_Hooke` - フックの法則に基づく剛性行列を計算。
  * `BBFE_elemmat_solid_mat_linear` - 線形要素の剛性行列を計算。
  * `BBFE_elemmat_solid_tensor_defgrad` - 変形勾配テンソルを計算。
  * `BBFE_elemmat_solid_mat_dispstr_tl` - Total Lagrangian法の変位-ひずみ行列を計算。
  * `BBFE_elemmat_solid_tensor_Green_Lagrange_mat_notation` - グリーン・ラグランジュテンソルを行列表現で計算。
  * `BBFE_elemmat_solid_tensor_second_Piora_Kirchhoff_mat_notation` - 二次ピオラ・キルヒホフテンソルを行列表現で計算。
  * `BBFE_elemmat_solid_mat_tl` - Total Lagrangian法の剛性行列を計算。
  * `BBFE_elemmat_solid_vec_inner_force_tl` - Total Lagrangian法の内部力ベクトルを計算。

***

### **FE\_manusol**

* **manusol.c / manusol.h**
  * `BBFE_manusol_calc_nodal_error_scalar` - 各節点でのスカラー量の誤差を計算。
  * `BBFE_manusol_set_bc_scalar` - スカラー量の境界条件を設定。

***

### **FE\_std**

* **integ.c / integ.h**
  * `BBFE_std_integ_calc` - 数値積分を実行。
  * `BBFE_std_integ_calc_volume` - 体積積分を実行。
  * `BBFE_std_integ_line_set_arbitrary_points` - 線積分の任意点を設定。
  * `BBFE_std_integ_rec_set_arbitrary_points` - 長方形領域での積分点を設定。
  * `BBFE_std_integ_tri_set_arbitrary_points` - 三角形領域での積分点を設定。
  * `BBFE_std_integ_hex_set_arbitrary_points` - 六面体領域での積分点を設定。
  * `BBFE_std_integ_tet_set_arbitrary_points` - 四面体領域での積分点を設定。
  * `BBFE_std_integ_tet_set_1points` - 四面体で1点積分を設定。
  * `BBFE_std_integ_tet_set_4points` - 四面体で4点積分を設定。
  * `BBFE_std_integ_tet_set_5points` - 四面体で5点積分を設定。

* **mapping.c / mapping.h**
  * `BBFE_std_mapping_void` - 空のマッピング処理。
  * `BBFE_std_mapping_calc_Jacobi_mat_3d` - 3Dのヤコビ行列を計算。
  * `BBFE_std_mapping_scalar` - スカラー量の空間マッピングを実行。
  * `BBFE_std_mapping_vector3d` - 3Dベクトル量の空間マッピングを実行。
  * `BBFE_std_mapping_scalar_grad` - スカラー量の勾配を計算。
  * `BBFE_std_mapping_vector3d_grad` - 3Dベクトル量の勾配を計算。
  * `BBFE_std_mapping_vector3d_div` - 3Dベクトル量の発散を計算。

* **shapefunc.c / shapefunc.h**
  * `BBFE_std_shapefunc_rec1st_get_val` - 1次の長方形形状関数の値を計算。
  * `BBFE_std_shapefunc_rec1st_get_derivative` - 1次の長方形形状関数の微分を計算。
  * `BBFE_std_shapefunc_tri1st_get_val` - 1次の三角形形状関数の値を計算。
  * `BBFE_std_shapefunc_tri1st_get_derivative` - 1次の三角形形状関数の微分を計算。
  * `BBFE_std_shapefunc_tri2nd_get_val` - 2次の三角形形状関数の値を計算。
  * `BBFE_std_shapefunc_tri2nd_get_derivative` - 2次の三角形形状関数の微分を計算。
  * `BBFE_std_shapefunc_hex1st_get_val` - 1次の六面体形状関数の値を計算。
  * `BBFE_std_shapefunc_hex1st_get_derivative` - 1次の六面体形状関数の微分を計算。
  * `BBFE_std_shapefunc_hex1st_get_surface` - 六面体形状関数の表面を取得。
  * `BBFE_std_shapefunc_tet1st_get_val` - 1次の四面体形状関数の値を計算。
  * `BBFE_std_shapefunc_tet1st_get_derivative` - 1次の四面体形状関数の微分を計算。
  * `BBFE_std_shapefunc_tet1st_get_surface` - 四面体形状関数の表面を取得。
  * `BBFE_std_shapefunc_tet2nd_get_val` - 2次の四面体形状関数の値を計算。
  * `BBFE_std_shapefunc_tet2nd_get_derivative` - 2次の四面体形状関数の微分を計算。
  * `BBFE_std_shapefunc_tet2nd_get_surface` - 四面体形状関数の表面を取得。

* **surface.c / surface.h**
  * `BBFE_std_surface_get_num_surfs_in_elem` - 要素内の表面数を取得。
  * `BBFE_std_surface_get_num_nodes_on_surf` - 表面上の節点数を取得。
  * `BBFE_std_surface_get_surface_node_3d` - 3Dの表面節点を取得。
  * `BBFE_std_surface_get_surface` - 表面情報を取得。
  * `BBFE_std_surface_hex1st_get_surface_node` - 1次六面体の表面節点を取得。
  * `BBFE_std_surface_tet1st_get_surface_node` - 1次四面体の表面節点を取得。
  * `BBFE_std_surface_tet2nd_get_surface_node` - 2次四面体の表面節点を取得。
  * `BBFE_std_surface_tet1st_get_surface` - 1次四面体の表面情報を取得。
  * `BBFE_std_surface_hex1st_get_surface` - 1次六面体の表面情報を取得。

***

### **FE\_sys**

* **memory.c / memory.h**
  * `BBFE_sys_memory_allocation_integ` - 数値積分用のメモリを確保。
  * `BBFE_sys_memory_free_integ` - 数値積分用メモリを解放。
  * `BBFE_sys_memory_allocation_shapefunc` - 形状関数用のメモリを確保。
  * `BBFE_sys_memory_free_shapefunc` - 形状関数用メモリを解放。
  * `BBFE_sys_memory_allocation_node` - 節点データ用のメモリを確保。
  * `BBFE_sys_memory_free_node` - 節点データ用メモリを解放。
  * `BBFE_sys_memory_allocation_elem` - 要素データ用のメモリを確保。
  * `BBFE_sys_memory_free_elem` - 要素データ用メモリを解放。
  * `BBFE_sys_memory_allocation_Dirichlet_bc` - Dirichlet境界条件用のメモリを確保。
  * `BBFE_sys_memory_free_Dirichlet_bc` - Dirichlet境界条件用メモリを解放。

* **read.c / read.h**
  * `BBFE_sys_read_node` - 節点データをファイルから読み込む。
  * `BBFE_sys_read_elem` - 要素データをファイルから読み込む。
  * `BBFE_sys_read_Dirichlet_bc` - Dirichlet境界条件をファイルから読み込む。

* **write.c / write.h**
  * `BBFE_sys_write_fopen` - 書き込み用ファイルを開く。
  * `BBFE_sys_write_vtk_shape` - VTK形式で形状データを出力。
  * `BBFE_sys_write_vtk_shape_with_disp` - 変位を含むVTK形式の形状データを出力。

***

### **libBB**

* **calc.c / calc.h**
  * `BB_calc_vec3d_dot` - 3Dベクトルの内積を計算。
  * `BB_calc_vec3d_cross` - 3Dベクトルの外積を計算。
  * `BB_calc_vec3d_length` - 3Dベクトルの長さを計算。
  * `BB_calc_vec3d_distance` - 2点間の距離を計算。
  * `BB_calc_vec3d_normal_vec` - ベクトルの法線ベクトルを計算。
  * `BB_calc_vec3d_copy` - ベクトルをコピー。
  * `BB_calc_mat3d_determinant` - 3D行列の行列式を計算。
  * `BB_calc_mat3d_inverse` - 3D行列の逆行列を計算。
  * `BB_calc_mat3d_copy` - 行列をコピー。

* **std.c / std.h**
  * `BB_std_calloc_1d_double` - 1次元のdouble型配列を確保。
  * `BB_std_free_1d_double` - 1次元のdouble型配列を解放。
  * `BB_std_calloc_2d_double` - 2次元のdouble型配列を確保。
  * `BB_std_free_2d_double` - 2次元のdouble型配列を解放。
  * `BB_std_calloc_3d_double` - 3次元のdouble型配列を確保。
  * `BB_std_free_3d_double` - 3次元のdouble型配列を解放。
  * `BB_std_read_file_get_val_double` - ファイルからdouble型の値を読み込む。
  * `BB_std_read_file_get_val_int` - ファイルからint型の値を読み込む。

* **vtk.c / vtk.h**
  * `BB_vtk_write_header` - VTK形式のヘッダーを書き込む。
  * `BB_vtk_write_points_3d` - 3Dの節点データをVTK形式で出力。
  * `BB_vtk_write_points_3d_with_disp` - 変位を含む3D節点データをVTK形式で出力。
  * `BB_vtk_write_cells` - 要素データをVTK形式で出力。
  * `BB_vtk_write_cell_types` - 要素タイプをVTK形式で出力。
  * `BB_vtk_write_point_vals_scalar` - スカラー量の節点値をVTK形式で出力。
  * `BB_vtk_write_point_vals_vector` - ベクトル量の節点値をVTK形式で出力。

***

## メッシュと境界条件の生成

BBFEフレームワークでは、解析対象のメッシュ生成や境界条件設定を、`util`ツールを用いて効率的に行います。主な利用方法を以下に説明します。

#### 1. メッシュ生成

メッシュ生成は、`meshgen`シリーズのツールを使って行います。以下に主要なメッシュ生成コマンドを紹介します。

* **六面体メッシュ生成 (meshgen\_hex)**
  ```bash
  ./meshgen_hex [x方向の要素分割数] [y方向の要素分割数] [z方向の要素分割数] [x方向の長さ] [y方向の長さ] [z方向の長さ]
  ```
  例えば、`meshgen_hex 10 10 10 1.0 1.0 1.0` と実行することで、10×10×10分割の各辺が1.0の六面体メッシュを生成し、節点データ`node.dat`と要素データ`elem.dat`が出力されます。

- **四面体メッシュ生成 (meshgen\_tet)**
  ```bash
  ./meshgen_tet [x方向の要素分割数] [y方向の要素分割数] [z方向の要素分割数] [x方向の長さ] [y方向の長さ] [z方向の長さ]
  ```

六面体メッシュを基に四面体要素を生成することが可能です。四面体要素は各六面体を6分割して得られます。

#### 2. 表面メッシュ抽出

生成されたメッシュの表面部分を抽出するために、`surf_conn`などのツールを使用します。

* **表面抽出 (surf\_conn)**
  ```bash
  ./surf_conn [ブロック長]
  ```

`surf_conn 3` のようにコマンドを実行すると、三次元のベクトル量に対する表面メッシュの情報（`surf.dat`）が出力されます。抽出された表面メッシュは、後に境界条件設定に利用されます。

#### 3. 表面の切り出し・除外

特定の領域から表面要素を抽出したり除外したりするには、`mesh_surf_extract`や`mesh_surf_remove`などのコマンドを使用します。

* **表面の切り出し (mesh\_surf\_extract)**

  ```bash
  ./mesh_surf_extract [x_min] [y_min] [z_min] [x_max] [y_max] [z_max]
  ```

指定した範囲の表面要素を抽出し、`surf.dat` として保存します。この操作により、解析対象の表面だけを取り出して境界条件を付加することが容易になります。

* **表面の除外 (mesh\_surf\_remove)**

  ```bash
  ./mesh_surf_remove [x_min] [y_min] [z_min] [x_max] [y_max] [z_max]
  ```

指定した範囲の表面要素を除外することで、不要な領域を取り除きます。抽出したい特定の部分を明確にする際に役立ちます。

#### 4. 境界条件の設定

メッシュに対して境界条件を設定するためには、`surf_dbc`や`surf_nbc`を使用します。

* **Dirichlet境界条件の付加 (surf\_dbc)**
  ```bash
  ./surf_dbc [ブロック長: n] [境界条件値 1] [境界条件値 2] ... [境界条件値 n]
  ```

表面の節点に対して指定した値のDirichlet境界条件を付加します。例えば、`surf_dbc 3 0.0 0.0 0.0` を実行することで、3次元ベクトル量に対する境界条件（{0.0, 0.0, 0.0}）を設定します。

* **Neumann境界条件の付加 (surf\_nbc)**

  ```bash
  ./surf_nbc [ブロック長: n] [境界条件値 1] [境界条件値 2] ... [境界条件値 n]
  ```

入力した表面上の節点に対して、指定した値のNeumann境界条件を付加します。このコマンドを用いることで、外力や熱流束などの境界条件を定義することができます。

* **境界条件ファイルのマージ (surf\_bc\_merge)**

  ```bash
  ./surf_bc_merge [境界条件ファイル 1] [境界条件ファイル 2]
  ```

  既存の2つの境界条件ファイルをマージして、全体の境界条件ファイルを作成します。例えば、異なる部分に設定した速度と圧力の条件を1つのファイルにまとめる場合に使用します。

#### 5. 使用例

例えば、`meshgen_hex`を使用してメッシュを生成し、`surf_conn`で表面メッシュを抽出し、`surf_dbc`を用いてDirichlet境界条件を付加する一連の流れを実行することで、解析の準備を整えることができます。

```shellscript
#!/bin/bash

# メッシュ生成
meshgen_hex 10 10 10 1.0 1.0 1.0

# 表面メッシュの抽出
surf_conn 3

# 表面の一部にDirichlet境界条件を付加
surf_dbc 3 0.0 0.0 0.0
```

このスクリプトにより、10×10×10のメッシュ生成、表面メッシュの抽出、Dirichlet境界条件の付加が一気に実行されます。

## Makefileによるコンパイル

BBFEフレームワークの使用には、解析ソルバーやライブラリのコンパイルが必要です。この過程ではMakefileを利用します。Makefileを使うことで、必要なコンパイルオプションやライブラリのパス設定を容易に管理できます。

#### 1. Makefileの概要

BBFEフレームワークのソルバーをコンパイルするためのMakefileは、以下のような役割を持ちます。

* 必要なライブラリのインクルードディレクトリやリンクディレクトリの指定
* ソルバーのオブジェクトファイルを生成し、それをリンクして実行ファイルを生成する

Makefileの主な構成要素には以下があります：

* **INCLUDES**: ヘッダファイルの場所を指定します。例えば、BBFEライブラリや`monolis`に関連するヘッダファイルのパスを指定します。
* **LIBS**: リンクするライブラリのパスを指定します。
* **TARGET**: 生成する実行ファイルの名前を指定します。

例えば、`fluid_fs`ソルバーのMakefileは次のような内容になります:

```Makefile
CC = mpic++
FC = mpif90

CFLAGS = -g -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../submodule/monolis/include
LIBS = -L/usr/local/lib -lm -lstdc++ -L../../submodule/monolis/lib -lmonolis -lmetis

INCLUDES_BB = -I../../include
LIBS_BB = -L../../lib -lBB -lBBFE_std -lBBFE_sys -lBBFE_elemmat -lBBFE_manusol

TARGET1 = fluid_fs
OBJS1 = fluid_fs.o fluid_core.o

.SUFFIXES: .c .cpp .o

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(FC) $(CFLAGS) $(INCLUDES) $(INCLUDES_BB) $(LDFLAGS) -o $@ $^ $(LIBS) $(LIBS_BB)

clean:
		-rm -f $(OBJS1) $(TARGET1)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES_BB) $(LIBS_BB)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES_BB)  $(LIBS_BB)
```

#### 2. Makefileの変更点

BBFEフレームワークのライブラリを新しい場所に移動した場合や、ソルバーを別の場所に移植した場合には、Makefileの中で指定されたパスを適切に修正する必要があります。例えば、以下のように設定を変更します​。

* INCLUDESの変更: `-I../../submodule/monolis/include`などのパスを、新しいライブラリの場所に合わせて修正します。
* LIBSの変更: `-L../../submodule/monolis/lib` などのパスも同様に変更が必要です。
* BBFEライブラリ関連の設定: `INCLUDES_BB` や `LIBS_BB` についても、新しいパスに応じて変更する必要があります。

#### 3. コンパイルの手順

Makefileを使ってソルバーをコンパイルする際の基本手順は以下の通りです。

1. ライブラリのインストール:

   * 最初にBBFEライブラリをインストールします。BBFE関連のシェルスクリプト（例：`install.sh`）を実行することで必要なライブラリが`lib`ディレクトリに生成されます。

2. ソルバーのコンパイル:

   * コンパイルしたいソルバーのディレクトリに移動し、`make`コマンドを実行します。
   * これにより、Makefileに記述された内容に基づいてソルバーがコンパイルされ、実行ファイルが生成されます。

#### 4. `solve_mat`の使用とMakefile

`monolis`が使いにくい場合、代替として`solve_mat`を使用することができます。`solve_mat`を使用するには、`Makefile`内での適切な設定が必要です。`solve_mat`は`monolis`とは異なる疎行列格納形式を用いているため、それに応じたパスやライブラリの設定を行う必要があります。
