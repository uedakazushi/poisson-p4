# AB原点四次のcohomology上界の独立監査

2026年11月2日。[四次主証明](AB_ORIGIN_DEGREE_FOUR_20261101.md)の
H¹(E4)=0、dim D4=38570を独立に検証した。
評価minorによる座標環B4の下界は別担当の監査であり、このレビューでは再実装しない。

## 1. Gröbner基底を使わないLevi character

VをLeviの標準二次元表現、E=Sym²Vとすると、
coreの線形表示はEnd(E)⊗det Vのtrace-zero部分である。
rank-one行列coneの次数n部分のcharacterは

    (det V)^n Sym^n(E) ⊗ Sym^n(E*)

である。traceは重みdet Vの非零因子なので、Yの次数nのcharacterはこれから

    (det V)^n Sym^(n−1)(E) ⊗ Sym^(n−1)(E*)

を引いて得られる。n=0では後項を零とする。
これを独立に列挙し、dim Y_n=(n+1)^3と全weight重複を得た。
主producerのcore Gröbner基底・先頭単項式やbackendはimportしない。

verticalの重み(3,1),(2,2),(1,3)の対称冪を掛け、各Levi既約型へ分解した。
ρ=(2,1)を加えた(x,y)に対して、x−y,y,x+y,xの負符号数からBott次数を計算する。
符号・双対bundle規約は[三次の独立レビュー](AB_ORIGIN_DEGREE_THREE_REVIEW_20261101.md)
で一次資料Fonarev Theorem3.2を確認したものと同じである。
四次の全Bott表が主reportと一致し、H⁰=41621、H¹=3051、H²以上は零となった。

## 2. P-cyclic vectorの評価制約

W_(a,b)のlowest vector y^(a−b)(y∧w)^bの射影G軌道は、b>0ならLagrangian旗、
b=0なら直線の旗である。固定q平面に横断する旗はP軌道の開集合をなす。
従ってこのvectorのP線形包はG線形包と等しく、既約表現全体を生成する。
この議論はb=2の場合にも適用され、b=1だけの計算を外挿していない。

G-mapの固定平面での評価はinjectiveにP-mapへ写る。
評価が零ならGの推移性により全点で零だからである。
P-mapはcyclic vectorの像だけで定まるので、そのweight空間における必要な
annihilator条件からH⁰内の重複数を抑えられる。

今回はpush/extract helperを使わず、二つの微分作用を直接導いた。
Levi upperはS(x,z),H(x,z),g(x,z)にx→x−uzを作用させる。
radical w→w+uzはS,Hを固定し、

    δg=H_xx−(z/4)S_x

を与える。従って(δg0,δg1,δg2)=(12h0,6h1−s0/2,2h2−s1/4)である。
lowest vectorはLevi upperに固定され、radical作用はu次数高々bなので、
像には前者の一階微分零と後者の(b+1)階微分零が必要である。

trace-zero rank-oneの稠密chart N=(1,a,b)^t(−ad−be,d,e)へ代入し、
weight空間の独立性と両作用の制約行列を有理数で再構成した。

| 型 | weight空間次元 | 制約階数 | Hom上界 |
|---|---:|---:|---:|
| (12,2) | 1 | 1 | 0 |
| (11,1) | 3 | 2 | 1 |
| (10,0) | 7 | 5 | 2 |
| (12,0) | 1 | 1 | 0 |

Yの関係を忘れた自由多項式環での独立性や、必要条件からP-mapの存在を推論していない。

## 3. H¹消滅と範囲

五段filtrationの各graded pieceにH²以上がないため、E4にもH²以上はない。
H¹に現れ得る四つの型について、独立に得たH⁰重複数上界はEuler重複数に等しい。
H⁰−H¹=Eulerと非負性により全H¹重複数が零となる。
従ってH¹(E4)=0、dim D4=41621−3051=38570である。

このレビューだけで全次数の生成性や原点の正規性を結論しない。
それらには別の有限生成基準が必要である。

```
.venv/bin/python code/ab_origin_degree_four_cohomology_review_20261101.py --output /tmp/ab-fourth-cohomology-review.json
```
