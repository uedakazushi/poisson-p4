# S61 の成分対応と重み1部分の独立点検

2026-10-30。対象は `legacy/v02/p4_poisson_completeness_v2.tex` §3 と
同版 `code/build_symplectic_slice.py`、`code/reduce_symplectic_bianchi.py`。
このノートは現行研究用の補足であり、公表原稿・先行研究上の新規性の主張ではない。

## 1. 成分対応を使うための仮定

基礎体は C、対象は非零 Poisson bivector の射影化、PGL5 商を取る前の空間とする。
固定点 p=[1:0:0:0:0] における評価が非退化な開部分を U とする。
一般階数4の各既約成分は PGL5 に保たれる。実際、連結代数群は有限個の既約成分を
置換しない。各構造の非退化点を p に移せるので、各一般階数4成分は U と交わる。
したがって U の既約成分は求める一般階数4成分と一対一に対応する。

評価の値 [J] の安定化群は CSp4 で連結である。ここで
J=partial_1 wedge partial_2+partial_3 wedge partial_4 とする。
評価の像は GL4/CSp4、次元は16−11=5。U はこの等質空間上の associated bundle
であり、連結な安定化群は fibre の既約成分を置換しない。この議論では評価が非退化
であること、射影化していることの両方が必要である。定数項 J に一意に代表を選ぶと、
fibre は affine 係数次数1,2,3の24+60+36=120変数で表される。

4次元の perspective 群を除く際には、単なる接空間の商ではなく、有限の多項式的
正規化が必要である。下記の明示式がこの部分を直接説明する。正規化後に55変数を
重み順で消去する v2 の手続きは、非零の有理定数のみを分母とする。従って形式近傍に
限定されず代数的同型となる。残る154式から134式への変更は20個の多項式恒等式に
よる ideal の同一性であり、根基化ではない。

以上から、これらの exact elimination certificates が成立することを条件として、
S61 の成分対応と次元差5+4=9に論理上の問題は見つからなかった。
この読解は全134式の再生成の代替ではない。全式の実行検証は別途記録する。
v2 本文末の「既知40個だけ」を示すという歴史的記述は、現行の57の下界と整合するよう
「すべての極小成分を尽くす」と読む必要がある。S61 に正確に40成分しかないとはしない。

## 2. Perspective 変換の有限公式

以下は actual pushforward の符号規約で記す。v2 の
`linear_action_derivative` は Lie derivative を計算するので、無限小の符号は逆である。
ベクトル y=(y1,...,y4)^t、Euler 場 E=y、反対称係数行列としての J を用いる。
P=J+L+Q+C とし、各項の係数次数は0,1,2,3とする。射影延長条件により
C=E wedge B、B は二次ベクトル場である。a を4成分の列ベクトルとして

    g_a(y)=y/(1+a^t y),   f=a^t y,   h=1-f

と置く。この変換は x0'=x0+sum a_i x_i, x_i'=x_i に由来する。
g_a^{-1}(y)=y/h における微分は h(I−y a^t) である。任意の反対称行列 M に対して

    (I−y a^t) M (I−a y^t)=M+E wedge (M a)

が成立する。さらに M=C の場合、右辺は h C である。したがって

    (g_a)_*P = h^2 (J+E wedge J a)
                + h (L+E wedge L a)
                + Q+E wedge Q a+C.

特に次数別には

    L' = L−2fJ+E wedge J a,
    Q' = Q−fL+E wedge L a+f^2 J−2f E wedge J a,
    C' = C+E wedge Q a−f E wedge L a+f^2 E wedge J a.

これは有限変換そのものの公式で、一次近似ではない。係数は a と元の係数の多項式であり、
a に重み1を与えると各項の重みを保つ。g_a g_b=g_{a+b} で逆変換は g_{−a}。
元の各 Poisson 構造を保つ対称群をここで分類しているわけではない。

## 3. 16方向の自然な正規化

一次 bivector L に対し

    tau(L)=L_12+L_34=(1/2) tr(J^{-1} L)

を線形形式とする。v=J a なら (E wedge v)_12+(E wedge v)_34=−a^t y。
従って上の式から

    tau(L')=tau(L)−5a^t y.

よって a=(1/5) grad(tau(L)) により、**tau(L')=0 に一意に正規化できる**。
この断面は任意の線形 jet 上で定義され、局所化を必要としない。
Jacobi の次数0の条件 [J,L]=0 は24変数に対する独立な4線形条件であるから、
閉じた一次 bivector の空間 Z1 は20次元。perspective 方向

    T_a=−2(a^t y)J+E wedge J a

は [J,T_a]=0、tau(T_a)=−5a^t y を満たす。従って

    Z1 = {L : [J,L]=0, tau(L)=0} direct_sum {T_a : a in C^4},
    dim {L : [J,L]=0, tau(L)=0}=16.

具体的な射影は L -> L+T_{grad(tau(L))/5} である。これは16方向の座標に依存しない
実現を与える。保存済み S61 の特定の4係数を零にする gauge と同一ではないが、両者の
間は上の有限変換で多項式的に移れる。重み2変数はこの変更で重み1変数の二次式を伴う
ため、保存済み134式に対し一次変数だけを機械的に取り替えてはいけない。

A J A^t=lambda J を満たす A in CSp4 の作用を、定数項を保つよう

    P(y) -> lambda^{-1} A P(A^{-1}y) A^t

と正規化すると tau(L) は tau(L)(A^{-1}y) に移る。従って tau=0 の gauge は
CSp4 で保たれ、16方向上の残存作用は明示的な線形作用になる。ここでは任意の
16方向が高次に積分できることも、全 S61 の極小素分解も主張しない。

補足として、16次元の ambient jet 空間全体を次元11の CSp4 の有限個の軌道だけで
覆うことは次元から不可能である。この次元比較だけでは S61 からの重み1射影の像について結論しない。
その後の独立な toric 族の計算により、像が非空開集合を含むことは証明された
（`TORIC_WEIGHT1_DOMINANCE_20261030.md`）。従って有限軌道での被覆は実際の jet 像でも不可能。
ただし像が全16次元空間に一致するとは示していない。

## 4. 独立検証と次の実行単位

`code/s61_weight1_audit.py` は旧実装を import せず、SymPy 1.14.0 上で QQ の exact
行列計算と多項式恒等式を検証する。17項目すべて成功。24係数の順序、Jacobi 行列、
tau 行列、perspective 行列を
`reports/20261030-s61-weight1/finite_identity_certificate.json` に保存した。
対象は rank4、rank8、rank16 の線形計算、trace 変換、冪等射影、有限変換の行列恒等式、
radial 三次項の相殺である。数値標本から幾何学的結論を推論していない。

再現は新しい出力ファイルを指定する。

```sh
.venv/bin/python code/s61_weight1_audit.py \
  --output reports/新しい検証名/finite_identity_certificate.json
```

次の有限目標は、保存済み61変数から primitive gauge への変換を明示し、重み1の
CSp4 共変量で層別を始めることである。同時に既知族の一般パラメータをこの gauge に
移し、重み1射影の像を求めれば、ambient jet 表現のうち実際に分類へ関与する部分を
特定できる。有限個のサンプルの検算を一般パラメータの証明の代わりにはしない。
