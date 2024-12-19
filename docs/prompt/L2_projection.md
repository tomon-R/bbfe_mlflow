L2プロジェクションは、任意の関数 ( u \in L^2(\Omega) ) を有限要素空間 ( V\_h \subset L^2(\Omega) ) に射影する方法です。このプロセスは、以下のような最小化問題として定式化されます：

### 数式による定義

有限要素空間 ( V\_h ) において、( u\_h \in V\_h ) を求める：
\[
J(u\_h) := \frac{1}{2} | u\_h - u |^2\_{L^2(\Omega)} \to \min
]
ここで、( | \cdot |\_{L^2(\Omega)} ) は ( L^2 ) ノルムを意味します。

### 解の性質

* ( J(\cdot) ) は凸関数であるため、唯一の解が存在します。
* フレシェ微分を用いて、最適性条件は以下のように表されます：
  \[
  (u\_h - u, v\_h) = 0 \quad \forall v\_h \in V\_h
  ]
  これは次のように簡約できます：
  \[
  (u\_h, v\_h) = (u, v\_h) \quad \forall v\_h \in V\_h
  ]

### 離散化

有限要素法を用いた近似では、離散系は以下の形を取ります：
\[
M x = b
]

#### 要素の定義

* 質量行列 ( M )：
  \[
  m\_{ij} = (\phi\_j, \phi\_i)
  ]
  ここで、( {\phi\_i} ) は有限要素空間 ( V\_h ) の基底です。
* 右辺ベクトル ( b )：
  \[
  b\_i = (u, \phi\_i)
  ]
* ( u\_h ) の表現：
  \[
  u\_h = \sum\_i x\_i \phi\_i
  ]
  ここで、( x ) は自由度のベクトルです。

### まとめ

L2プロジェクションは、元の関数 ( u ) を有限要素空間 ( V\_h ) に最も近い形で表現する手法であり、質量行列と右辺ベクトルを用いて線形代数の問題として解くことができます。

***

この式をフレシェ微分という言葉を使わずに導出する方法を説明します。

### 問題の復習

( u\_h \in V\_h ) を、以下の最小化問題の解とします：
\[
J(u\_h) = \frac{1}{2} | u\_h - u |^2\_{L^2(\Omega)} \to \min
]
ここで、( | v |^2\_{L^2(\Omega)} = \int\_\Omega v^2 , dx ) です。

有限要素空間 ( V\_h ) は有限次元空間であり、関数 ( u\_h ) をその空間内で最適に選ぶ問題です。

***

### 導出のステップ

#### 1. 目的関数の展開

目的関数 ( J(u\_h) ) を展開します：
\[
J(u\_h) = \frac{1}{2} \int\_\Omega (u\_h - u)^2 , dx
]

この積分を展開すると：
\[
J(u\_h) = \frac{1}{2} \int\_\Omega u\_h^2 , dx - \int\_\Omega u\_h u , dx + \frac{1}{2} \int\_\Omega u^2 , dx
]

#### 2. 最小化条件

最小化問題では、目的関数 ( J(u\_h) ) が変化しない方向（すべての摂動方向 ( v\_h \in V\_h ) に対して）において微小な変化をゼロにする必要があります。

つまり、任意の ( v\_h \in V\_h ) に対して、次が成り立つようにします：
\[
\frac{d}{d\epsilon} J(u\_h + \epsilon v\_h) \Big|\_{\epsilon=0} = 0
]

#### 3. 摂動を考慮した ( J(u\_h + \epsilon v\_h) )

( u\_h + \epsilon v\_h ) を ( J ) に代入します：
\[
J(u\_h + \epsilon v\_h) = \frac{1}{2} \int\_\Omega (u\_h + \epsilon v\_h - u)^2 , dx
]

この式を展開します：
\[
J(u\_h + \epsilon v\_h) = \frac{1}{2} \int\_\Omega \left\[ (u\_h - u)^2 + 2\epsilon (u\_h - u)v\_h + \epsilon^2 v\_h^2 \right] , dx
]

#### 4. 微小変化の寄与

(\epsilon) に関して微分を取ります：
\[
\frac{d}{d\epsilon} J(u\_h + \epsilon v\_h) = \int\_\Omega (u\_h - u)v\_h , dx + \epsilon \int\_\Omega v\_h^2 , dx
]

ここで、(\epsilon = 0) を代入すると：
\[
\frac{d}{d\epsilon} J(u\_h + \epsilon v\_h) \Big|*{\epsilon=0} = \int*\Omega (u\_h - u)v\_h , dx
]

#### 5. 最適性条件

最適な ( u\_h ) を得るためには、上式が任意の ( v\_h \in V\_h ) に対してゼロである必要があります：
\[
\int\_\Omega (u\_h - u)v\_h , dx = 0 \quad \forall v\_h \in V\_h
]

これを内積記号を使って簡潔に書くと：
\[
(u\_h - u, v\_h) = 0 \quad \forall v\_h \in V\_h
]

***

### 結論

この式は、有限要素空間 ( V\_h ) において ( u\_h ) が ( u ) に最も近い関数であることを保証します。これは、目的関数の微小変化を考え、その変化をゼロにする条件から導かれます。
