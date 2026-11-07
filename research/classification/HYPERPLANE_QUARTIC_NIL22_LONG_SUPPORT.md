# H+既約quartic：nil22・長translation全源の非支配

2026-11-07。`E_reduced_nonnormal` の全nilpotent表の型7、
`Z=y∂x+∂y+w∂z` を扱う。
`CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json` の
全90係数線形問題から得た20自由係数源 `n0,...,n19` を使い、全Jacobi係数を課す。

**固定cofactorの全Poisson源の次元は9以下、全射影飽和像の次元は24以下である。**
既約quartic条件を課す前の源全体について成り立つので、全退化・可約境界も含む。
従って未知の25–27次元優越成分の一般支持にはならない。
個々の点の成分帰属やE側の全被覆、新成分の存在はこの命題から主張しない。

最終算術証明書は `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT_20261107.json`、
reusable runner の mode は `--nil22-long-support`。
依存する有限次元証明書は同じprefixの `CORE`、`CRITICAL_CORES`、`LOW_RANK_CORE` である。

## 1. 非線形coreと二段の線形fiber

全20変数を次の三群に分ける。

    core = (n2,n8,n11,n12,n14,n17,n18,n19),
    a = (n1,n7,n10,n13,n16),
    b = (n0,n3,n4,n5,n6,n9,n15).

coreだけに依存する全Jacobi係数を課した基底をTとする。
全解は、まずT上の二つの斉次線形式 `A a=0`、次にaを固定した五つの
affine線形式 `B b=rhs(a)` を満たす。
元のJacobi係数表の0始まり添字はそれぞれ `(2,8)` と `(13,23,9,14,15)` である。
従ってこれは全解を含む二段fiberの上界であり、任意の追加整合条件は次元を増やさない。

`q=n14` と置く。最初の行列は

    A = [ -2n17, -n18, n18, n11+q-2n19, -n11+q
           2n18, 2n11-4q-2, -3n11+2n19+1, 2n12+2n8, -2n12-6n8 ].

次の定数恒等式により、全core点でrank A≥1である。

    A[1,1]+6A[0,4]−2A[0,3]−2A[1,2] = −4.

第二の行列は、上記bの順で

    B = [ 0,      0,      0,       -3q-1, 0,          0,                  0
          2n17,   0,      0,       0,     n18,       -n18,               n11-2q
          n18,    0,      0,       0,     n11-q-1,   (-3n11+2n19+1)/2, -n12-3n8
          2n19+1, 0,      -2q-1,   -2n13, 2n12,      2n8,                2n2
          n16,    -q-1,   -n13,    0,     n10,       n7,                 n1 ].

全rhsも最終証明書に保存した。

## 2. 第二fiberの一様階数

Bの第2,3行と列 `(n0,n6,n9,n15)` でできる2行4列blockをCとする。
C=0なら、線形forcingは正確に

    n17=n18=0, n11=2, q=1, n19=5/2, n12=−3n8

となる。forcing行列は定数の非零6次小行列式を持つ。
ところがこの置換でcoreのJacobi係数第25番は3となり、Tの点では不可能である。
従って全T上でC≠0となる。

Bの行 `(0,3,4)`、列 `(n5,n4,n3)` の三角小行列式は

    −(q+1)(2q+1)(3q+1).

Cの二行はこの三列で零である。従って
`q∉{−1,−1/2,−1/3}` では、三角pivotの階数3にCの非零成分が一つ加わり、rank B≥4となる。

三つのcritical値でも、残る二つの三角pivotはそれぞれ定数小行列式
`2,−1/4,2/9` を持つ。Cの非零成分を加えて、全critical点でrank B≥3である。
これは標本の階数でも一般点での階数でもなく、全core点・全aに対する下界である。

## 3. 全coreと全criticalcoreの次元証明

全core上界と三つのcriticalcore、さらに
`q=−1/2, rank A≤1` の閉層について、QQ上の有限多項式証明書を作った。
最後の閉層ではAの全10個の2次小行列式を課しており、選んだ一つのminorだけへ制限していない。

各証明書は元のcore係数から得た多項式Gと、各Gを入力理想の生成元の組合せとして表すliftを保存する。
全lift恒等式を別にSymPy QQ多項式演算で検査した。
全Gのleading monomialが生成するmonomial idealの高さを、有限vertex cover再帰で計算した。
Gが入力理想に属すればこのmonomial idealは入力のinitial idealに含まれるため、
**Gが完全なGroebner basisであることの追加仮定なしで**上界が従う。

必要なleading monomialの台をradical monomial idealとして書くと次の通りである。

| 全core条件 | leading support idealの生成元 | core次元上界 |
| --- | --- | --- |
| 全T | n11,n18,n2n12n19,n2n17,n8n12,n8n14,n8n17,n8n19,n12n14,n12n17,n14n17,n14n19 | 2 |
| q=−1 | n8,n11,n14,n17,n18,n19,n2n12 | 1 |
| q=−1/2 | n11,n14,n18,n2n12n19,n2n17,n8n12,n8n17,n8n19,n12n17 | 2 |
| q=−1/3 | n8,n11,n14,n17,n18,n19,n2n12 | 1 |
| q=−1/2, rank A≤1 | n2,n8,n11,n14,n17,n18,n19 | 1 |

例えば全Tの行ではn11,n18を除く六変数のどの三つにも、表の台が一つ含まれる。
従って自由に残せる変数は高々二つである。他の行も同じ有限組合せ論で確認できる。
反転するcore元はなく、全零・交差・特殊点を含む。

## 4. 全固定源の次元は9以下

上の全階数・全core上界を組み合わせる。

| 層 | core上界 | a-fiber上界 | b-fiber上界 | 全固定源上界 |
| --- | ---: | ---: | ---: | ---: |
| qが三critical値以外 | 2 | 4 | 3 | 9 |
| q=−1 | 1 | 4 | 4 | 9 |
| q=−1/3 | 1 | 4 | 4 | 9 |
| q=−1/2, rank A=2 | 2 | 3 | 4 | 9 |
| q=−1/2, rank A=1 | 1 | 4 | 4 | 9 |

rank A=0は第1節の定数恒等式で不可能であり、表は全点を覆う。
rhsの不整合でfiberが空になることは上界を強くするだけである。

## 5. G側の全非支配と未解決scope

固定線形部Jordan(2,2)の共役軌道は8次元であり、長translationの非零quotientは4次元の開集合である。
全affine作用 `(B,t)↦([B,A],Bc−At)` の20行20列行列を正確に保存し、階数12も計算した。
超平面の選択4次元を加え、射影scalarを除くと

    dim(saturated source) ≤ 9+12+4−1 = 24.

この型の円錐は非零scalar倍を保つ。実際、
`diag(λ²,λ,λ,1)` はZをλZへ送る。
[優越成分の幾何](DOMINANT_COMPONENT_GEOMETRY.md)第1節で未知G成分は25,26,27次元なので、
この全型はその一般支持にならない。
全源のどの点がどの既知成分に入るかというE側帰属は未解決として残す。

計算には公開環境の既設Singularと共有libraryを再利用した。
最初の全20変数liftstdは80秒で上限となり、証明を得なかった。
`...NIL22_LONG_RESOURCE_LIMIT_20261107.json` はその失敗と全未解決入力を保存している。
確定証明はこの失敗に依存せず、全core・全critical層・一様二段fiberへ数学的に経路を変えて得た。
