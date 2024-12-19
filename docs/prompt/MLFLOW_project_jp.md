# BBFE MLFLOW PROJECT

BBFE MLFLOW PROJECTは、C言語用有限要素法解析フレームワーク（BBFE）を用いて開発された、二層流れ解析ソルバーです。本プロジェクトでは、Navier-Stokes方程式を基にした流体解析に加え、レベルセット法を用いた自由表面追跡を実現しています。

***

## 目次

- [BBFE MLFLOW PROJECT](#bbfe-mlflow-project)
  - [目次](#目次)
  - [解析手法](#解析手法)
    - [非圧縮性粘性流れの安定化有限要素法](#非圧縮性粘性流れの安定化有限要素法)
    - [レベルセット法を用いた自由表面の表現](#レベルセット法を用いた自由表面の表現)
  - [ソルバー構成](#ソルバー構成)
    - [Fractional Step法（mlflow\_fs）](#fractional-step法mlflow_fs)
    - [SUPG法を用いた安定化有限要素法（mlflow\_fs\_sups）](#supg法を用いた安定化有限要素法mlflow_fs_sups)
  - [解析事例](#解析事例)
    - [ダムブレイク問題](#ダムブレイク問題)
    - [気泡上昇流れ問題](#気泡上昇流れ問題)
    - [スロッシング解析](#スロッシング解析)
  - [Valuesに定義する変数](#valuesに定義する変数)
  - [プロジェクト特有のライブラリ群](#プロジェクト特有のライブラリ群)
    - [util](#util)
      - [`levelset.dat`ファイルの作成手順](#levelsetdatファイルの作成手順)
        - [**1. コマンドを使用して生成**](#1-コマンドを使用して生成)
        - [**2. プログラムの主な動作**](#2-プログラムの主な動作)
        - [**3. `levelset_gen` の利用例**](#3-levelset_gen-の利用例)
        - [**4. 注意点**](#4-注意点)
    - [関数群](#関数群)
      - [**BBFE\_fluid**](#bbfe_fluid)
      - [**BBFE\_elemmat**](#bbfe_elemmat)
      - [**BBFE\_mlflow**](#bbfe_mlflow)
      - [**mlflow\_utils**](#mlflow_utils)

***

## 解析手法

### 非圧縮性粘性流れの安定化有限要素法

流体は非圧縮性の粘性流体であり，Navier-Stokes方程式と連続の式が用いられる。fractional step法を採用した分離型解法ソルバ(mlflow\_fs)とGLS法を採用した直接型解法ソルバ(mlflow\_fs\_sups)がある。

### レベルセット法を用いた自由表面の表現

* **アルゴリズム**: レベルセット関数を用いて界面を追跡し、Navier-Stokes方程式を解く。レベルセット関数は界面の移流方程式に基づき更新。再初期化によって関数の性質を維持。体積保存性を確保する補正手法を採用。

* **表面張力の計算**: CSFモデルを用い、近似デルタ関数による離散化で界面付近の張力を計算。

* **再初期化手法**：ハミルトン・ヤコビ方程式を用いてレベルセット関数を再初期化。  特異点の処理や収束性の問題に対処する安定化項を導入。

* **体積補正法**：気液間の体積保存のため、界面体積を計算し、修正量を加えて体積誤差を補正。

## ソルバー構成

### Fractional Step法（mlflow\_fs）

* 分離型解法を採用。
* 速度場と圧力場を段階的に解くことで計算効率を向上。

### SUPG法を用いた安定化有限要素法（mlflow\_fs\_sups）

* 直接型解法を採用。
* 非圧縮性を保つ高精度な安定化手法。

## 解析事例

### ダムブレイク問題

水と空気の密度差を考慮したベンチマーク問題。
衝撃捕捉項の有無による安定性の検証。

### 気泡上昇流れ問題

表面張力と粘性が支配的なケースに分けた解析。
解析結果を文献と比較し、レベルセット法の適用性を確認。

### スロッシング解析

矩形容器内の自由表面の動きをALE法を用いて解析。
衝撃捕捉項の有効性を確認。

## Valuesに定義する変数

`Values`構造体は、MLFLOWソルバーの計算条件や物理量を管理するために使用されます。以下は、`Values`構造体に定義される主な変数です。

```c
typedef struct {
    int    num_ip_each_axis;       // 各軸方向の積分点の数
    double mat_epsilon;            // 行列計算の収束判定に用いる閾値
    int    mat_max_iter;           // 行列計算の最大反復回数

    double dt;                     // 時間ステップ幅
    double finish_time;            // 解析終了時間
    int    output_interval;        // 出力間隔（ステップ数）

    double   density_l, density_g; // 液相と気相の密度
    double   viscosity_l, viscosity_g; // 液相と気相の粘性係数
    double*  gravity;              // 重力加速度ベクトル
    double*  accel_amp;            // 外力加速度の振幅
    double*  accel_angle_vel;      // 外力加速度の角速度
    double*  accel_inertia;        // 慣性加速度
    double** v_mesh;               // メッシュ移動の速度ベクトル

    double volume_g_init;          // 気相の初期体積
    double volume_l_init;          // 液相の初期体積

    double** v;                    // 流体の速度ベクトル
    double*  p;                    // 圧力
    double*  density;              // ノードごとの密度
    double*  viscosity;            // ノードごとの粘性係数

    double* levelset;              // レベルセット関数（気相と液相の境界を表現）
    double* levelset_tmp;          // レベルセット関数の一時変数
    double* heaviside;             // ヘビサイド関数（レベルセットを基にした体積分率）

    double   surf_tension_coef;    // 表面張力係数
    double** surf_tension;         // 表面張力ベクトル
    double** grad_phi;             // レベルセット関数の勾配
    double   size_interface;       // インターフェースの厚み

    double dt_reinit;              // レベルセットの再初期化の時間ステップ幅
    double epsilon_reinit;         // 再初期化の安定化パラメータ ε
    double delta_reinit;           // 再初期化のデルタ関数パラメータ
    double alpha_reinit;           // 再初期化の制御パラメータ α
    int    max_iter_reinit;        // 再初期化の最大反復回数

    int output_option;             // 出力オプション
    int ale_option;                // ALE（Arbitrary Lagrangian-Eulerian）法の適用オプション

    int* measurement_node_id;      // 測定対象ノードのIDリスト
    int  measurement_num_nodes;    // 測定対象ノードの総数
} VALUES;

```

## プロジェクト特有のライブラリ群

MLFLOWソルバーは、BBFEフレームワークを基盤とし、特有のライブラリを拡張して利用します。以下はプロジェクトに特化したライブラリ群の概要です。

### util

#### `levelset.dat`ファイルの作成手順

`levelset.dat` ファイルは、解析領域の初期レベルセット関数を定義するためのファイルです。このファイルは、シミュレーションの初期条件を設定するために使用されます。以下に具体的な作成手順を示します。

##### **1. コマンドを使用して生成**

`levelset_gen` プログラムを使用して、以下のコマンドで生成します。

```bash
./levelset_gen SIMULATION_ID [x_min] [y_min] [z_min] [x_max] [y_max] [z_max]
```

* **SIMULATION\_ID**: シミュレーションタイプを指定。
  * `1`: ダム崩壊 (Dambreak)
  * `2`: 流体キャビティ (3D Cavity)
  * `3`: 気泡上昇 (Rising Bubble, 2D)
  * `4`: 気泡上昇 (Rising Bubble, 3D)
* **\[x\_min], \[y\_min], \[z\_min], \[x\_max], \[y\_max], \[z\_max]**:
  * 解析ドメインの境界を指定 (必要に応じて使用)。

例: ダム崩壊シミュレーションの初期条件を生成する場合

```bash
./levelset_gen 1 0.0 0.0 0.0 1.0 1.0 1.0
```

##### **2. プログラムの主な動作**

1. **引数の処理**
   * `args_manager_levelset` 関数がコマンドライン引数を解析し、入力/出力ディレクトリやファイル名を設定します。
   * SIMULATION\_ID に応じて、シミュレーションタイプと解析ドメインを設定します。

2. **距離関数の計算**
   * `calc_distance` 関数が、指定した境界条件に基づいて各メッシュ点からの距離を計算します。
   * 距離計算は、対象となるシミュレーション (ダム崩壊や気泡上昇) の幾何学的特性に基づいて実施。

3. **レベルセット関数の生成**
   * `write_levelset_file` 関数が距離情報を基に `levelset.dat` を生成します。
   * シミュレーションの初期条件に応じて、正の距離（液相）と負の距離（気相）を設定。

4. **出力**
   * 出力ファイル `levelset.dat` は以下の形式になります:
     ```
     #levelset
     <ノード数> 1
     <ノード1の値>
     <ノード2の値>
     ...
     ```

##### **3. `levelset_gen` の利用例**

* **ダム崩壊シミュレーション**
  ```bash
  ./levelset_gen 1 0.0 0.0 0.0 1.0 1.0 1.0
  ```
  境界 `(0,0,0)` から `(1,1,1)` までの範囲でレベルセットを計算。

* **2D 気泡シミュレーション**
  ```bash
  ./levelset_gen 3
  ```
  中心が `(0.5,0.5)`、半径 `0.25` の円形気泡を定義。

* **3D 気泡シミュレーション**
  ```bash
  ./levelset_gen 4
  ```
  中心が `(0.5,0.5,0.5)`、半径 `0.25` の球形気泡を定義。

##### **4. 注意点**

* 必要に応じて、入力ファイル (`node` ファイルや `element` ファイル) を準備してください。
* 実行前に、適切なディレクトリとファイル名を指定してください。
* 初期値の設定や計算範囲を正確に指定することで、解析結果の精度が向上します。

### 関数群

* **BBFE\_fluid**
  * 流体解析に特化したモジュール。
  * 主な機能：
    * 流体の速度 (`v`) と圧力 (`p`) の更新。
    * 非圧縮性 Navier-Stokes 方程式の安定化計算 (SUPG 法、GLS 法)。
    * メッシュ移動に伴う速度場の更新 (ALE 法)。
    * 質量保存を考慮した体積補正。
* **BBFE\_elemmat**
  * 有限要素法の基礎計算を提供するライブラリ。
  * 主な機能：
    * 要素剛性行列の生成。
    * ヤコビアン行列の計算。
    * 要素間の相互作用計算（衝撃捕捉項の追加を含む）。
* **BBFE\_mlflow**
  * BBFE MLFLOW プロジェクトの主要な拡張機能を提供するライブラリ。
  * 主な機能：
    * 界面追跡用の体積補正計算 (`volume_correction`)。
    * 表面張力計算（CSF モデルを使用）。
    * レベルセット関数とヘビサイド関数を利用した物性値更新。
* **mlflow\_utils**
  * MLFLOW の解析を補助するためのユーティリティ。
  * 主な機能：
    * 設定ファイルの読み込み (`cond.dat`, `levelset.dat`)。
    * 結果ファイルの書き出し (`vtk_export` と連携)。
    * 自由表面を表現するための補助関数。

#### **BBFE\_fluid**

* **流体解析に特化したモジュール**。有限要素法を使用した流体解析に必要な基本操作とSUPG法・GLS法のサポートを提供。

* **主な機能**:
  1. **データの初期化と設定**:
     * **`BBFE_fluid_pre`**:
       * 節点データ (`node.dat`) と要素データ (`elem.dat`) を読み込む。
       * 積分点の設定と形状関数の初期化を行う。
     * **`BBFE_fluid_set_basis`**:
       * 要素の形状（四面体または六面体）に応じて積分点や形状関数の導出。
       * 1次形状関数（1st-order）を用いて解析精度を確保。

  2. **速度・圧力の更新**:
     * **`BBFE_fluid_renew_velocity`**:
       * 計算結果ベクトルから速度場 (`v`) を抽出し、ノードごとに割り当てる。
     * **`BBFE_fluid_sups_renew_velocity`**:
       * SUPG法で計算された結果に基づき、速度場を更新。
     * **`BBFE_fluid_sups_renew_pressure`**:
       * 圧力 (`p`) を計算結果ベクトルから抽出し、密度 (`density`) を考慮してノードごとに設定。

  3. **メモリ管理**:
     * **`BBFE_fluid_finalize`**:
       * 積分点や形状関数、ノード・要素データなどの動的メモリを解放。

  4. **解析条件の設定**:
     * **`BBFE_fluid_get_directory_name`**:
       * コマンドライン引数からディレクトリ名を取得し、解析ファイルの管理を容易にする。

* **特徴的なポイント**:
  * **要素形状に応じた柔軟性**:
    * 四面体（tetrahedron）および六面体（hexahedron）に対応し、それぞれの形状に最適な積分点と形状関数を自動的に設定。
  * **Navier-Stokes 方程式のサポート**:
    * 非圧縮性流体解析において、SUPG法やGLS法の適用に必要な速度・圧力の計算を効率化。
  * **有限要素データ構造の統合**:
    * 節点データや要素データを統合的に管理し、解析の柔軟性を向上。

#### **BBFE\_elemmat**

* **有限要素法の基礎計算を提供するライブラリ**。解析の基盤となる要素剛性行列やベクトルの計算をサポート。

* **主な機能**:
  1. **要素剛性行列の生成**:
     * **`BBFE_elemmat_fluid_sups_mat`**:
       * SUPG法（Streamline Upwind Petrov-Galerkin）やPSPG法（Pressure Stabilizing Petrov-Galerkin）を用いて、流体の安定化された剛性行列を生成。
       * ALE（Arbitrary Lagrangian-Eulerian）法によるメッシュ速度も考慮。
     * **`BBFE_elemmat_mat_levelset`**:
       * レベルセット方程式のための要素行列を生成。
       * SUPG法を基にした安定化計算を含む。

  2. **ヤコビアン行列の計算**:
     * **`BBFE_elemmat_mat_CLSM_reinitialize`**:
       * 保存的レベルセット法（Conservative Levelset Method）の再初期化処理に必要なヤコビアン行列を計算。

  3. **要素間の相互作用計算**:
     * **`BBFE_elemmat_fluid_sups_vec`**:
       * 外力（重力、表面張力、慣性加速度）を考慮した流体の右辺ベクトルを生成。
       * 衝撃捕捉（Shock Capturing）や流れ場の安定化をサポート。
     * **`BBFE_elemmat_vec_surface_tension`**:
       * 表面張力を計算し、界面の張力ベクトルを生成。
     * **`BBFE_elemmat_vec_levelset_reinitialize`**:
       * レベルセット法の再初期化のためのベクトルを生成。

  4. **その他の計算**:
     * **`elemmat_supg_coef`**:
       * SUPG法に基づく安定化係数を計算。
     * **`BBFE_elemmat_vec_grad_phi_L2_projection`**:
       * レベルセットの勾配に基づくL2射影ベクトルを計算。
     * **`BBFE_elemmat_mlflow_shock_capturing_coef`**:
       * 衝撃捕捉係数（Shock Capturing Coefficient）を計算。

* **特徴的なポイント**:
  * **安定化手法の統合**:
    * SUPG法、PSPG法、LSIC（Least Squares Interpolation Constraint）法を組み合わせた高精度な数値計算。
  * **衝撃捕捉の考慮**:
    * 流体の衝撃波や急激な勾配に対応するため、衝撃捕捉項を要素剛性行列およびベクトルに追加。
  * **表面張力とレベルセットの統合**:
    * 表面張力と保存的レベルセット法の再初期化を効率的に処理する設計。

#### **BBFE\_mlflow**

* **BBFE MLFLOW プロジェクトの主要な拡張機能を提供するライブラリ**。

* **主な機能**:
  1. **界面追跡と体積補正**:
     * **`BBFE_mlflow_renew_vals_by_levelset`**:
       * レベルセット関数を使用して、流体内の物性値（密度や粘性など）を界面に沿って更新。
       * 各ノードにおいて、レベルセット関数の値に基づき、界面の局所物性値を線形補間。
     * **`BBFE_mlflow_convert_levelset2heaviside`**:
       * レベルセット関数をヘビサイド関数に変換して、界面のスムーズな遷移を表現。

  2. **表面張力計算**:
     * **`BBFE_mlflow_clear_surface_tension`**:
       * 表面張力ベクトルの初期化。
     * 他のライブラリと連携し、表面張力の効果を界面近傍に適用。

  3. **メッシュ操作と動的更新**:
     * **`BBFE_mlflow_renew_mesh_velocity`**:
       * 慣性加速度を考慮し、各ノードのメッシュ速度を更新。
     * **`BBFE_mlflow_renew_mesh_position`**:
       * メッシュの速度と時間ステップを基に、各ノードのメッシュ座標を更新。
     * **`BBFE_mlflow_renew_acceleration`**:
       * 加速度ベクトルを時間依存で更新。振動運動や外部慣性力をシミュレート可能。

  4. **レベルセット関数操作**:
     * **`BBFE_mlflow_renew_levelset`**:
       * レベルセット関数の値を時間発展させて更新。
     * **`BBFE_mlflow_convert_levelset2CLSM`**:
       * 保存的レベルセット法（Conservative Levelset Method, CLSM）用にレベルセット関数を変換。
       * 平滑な界面遷移を表現するために指数関数変換や双曲線正接関数を使用。

* **特徴的なポイント**:
  * **界面追跡**:
    * ヘビサイド関数や保存的レベルセット法を用いて、界面を精密にトラッキング。
    * 界面近傍での物性値遷移をスムーズに計算。
  * **メッシュ移動**:
    * 慣性効果や加速度を考慮して、動的なメッシュ更新を実現。
  * **表面張力の統合**:
    * 他のBBFEモジュールと統合して、界面物理現象（特に表面張力）の計算を効率化。

#### **mlflow\_utils**

* **MLFLOW の解析を補助するためのユーティリティライブラリ**。

* **主な機能**:
  1. **設定ファイルの読み込みと結果ファイルの出力**:
     * 設定ファイル (`cond.dat`, `levelset.dat` など) を解析で利用するために読み込む。
     * 結果ファイルの出力:
       * **`output_result_dambreak_data`**:
         * ダム崩壊シミュレーションの解析結果を出力。
         * 水位や自由表面の位置を CSV 形式で保存。
       * **`output_result_bubble_data`**:
         * 気泡上昇シミュレーションの解析結果を出力。
         * 気泡の平均高さ、速度、球形度、サイズを保存。
       * **`output_result_sloshing_data`**:
         * 液面振動（スロッシング）の解析結果を出力。
         * 自由表面の位置を記録。

  2. **自由表面に関連する補助機能**:
     * **自由表面の位置計算**:
       * 結果ファイル出力時に自由表面のゼロ交差を計算。
       * 各測定点でのレベルセット関数の値を利用して表面位置を算出。
     * **測定点の管理**:
       * **`count_mlflow_measurement_node`**:
         * 測定点（例えば特定の座標点における自由表面）の数を計算。
       * **`set_mlflow_measurement_node`**:
         * 測定点の ID を設定。

  3. **統計的計算とデータ変換**:
     * **`calc_data_bubble`**:
       * 気泡上昇における速度や高さ、体積の平均値を計算。
     * 並列環境での MPI 通信をサポートし、全体の統計を効率的に計算。
     * 計算結果は物理モデルやメッシュ構造に基づき精密に評価される。

* **特徴的なポイント**:
  * **モジュール全体の補助**:
    * MLFLOW のシミュレーション結果を効率的に処理し、後処理や可視化のためのデータを生成。
  * **MPI 対応**:
    * 並列計算環境におけるデータ収集と統計処理をサポート。
  * **結果データの汎用性**:
    * CSV ファイル形式で出力されるため、他の解析ツールや可視化ソフトと簡単に連携可能。
