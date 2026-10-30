# 一般トーリック族とalmost-abelian族のS61座標

2026-10-30。原資料v02の模型に対する現行の座標表示。文献上の新規性は主張しない。

## トーリック族の全6パラメータを保つchart

4×4反対称行列Lambdaの上三角成分を順に(a,b,c,d,e,f)とし、
Delta=af−be+cdと置く。開条件はa Delta≠0。
affine bracket `{z_i,z_j}=lambda_ij z_i z_j` は全パラメータでPoissonである。
可換なベクトル場 `z_i partial_i` の外積の定数係数線形結合だからである。
非退化条件Delta≠0の下で、実際の射影退化因子は5座標超平面の和となる。
この族は旧カタログのT01（無smoothing edge）に属する。
ここでLambdaはPoisson係数で、旧CSVにあるbiresidue行列そのものではない。

点z=(1,1,1,1)を原点に移し、接空間の新座標y=A(z−1)を取る。ただし

```text
A = [ 1       0        0       0       ]
    [ 0       1/a      0       0       ]
    [ d/a    -b/a      1       0       ]
    [ e/Delta -c/Delta 0       a/Delta ]
```

直接計算でA Lambda A^t=J、det(A)=1/Delta。したがって定数項がJとなる。
斉次二次双ベクトルの実際のpushforwardは `M P(M^-1 x) M^t` で行う。
続いてv02の4個のperspective座標を消す変換の係数は
`(1/3,a/3,0,0)` となる。Euler垂直成分を除き、保存済み自由係数を読むと、
61個の有理関数が得られる。

一般パラメータでの61座標は
`reports/S61_FAMILIES_20261030/toric_T01_map.json`、生成・検証は
`code/verify_s61_families.py`。同スクリプトは55復元式との一致、134式、
復元後のprojective Jacobiを有理関数体QQ(a,b,c,d,e,f)上で検証する。
分母を払った比較は、指定開集合上の恒等式を証明する。

この一つの写像は6パラメータであり、S61内のT01成分全体の次元16と同じではない。
残存するSp4作用を加えた一次部分への写像は、別の非零minorによる優勢性証明で扱う。
この点を無視して6変数だけで全成分を尽くしたとしない。

除外したa=0, Delta≠0は他の非零成分をpivotとするsymplectic chartで扱えるが、
その全chart atlasは今回実装していない。Delta=0はこの非退化評価chartの対象外。

## Almost-abelian族

v02の既存族

```text
P = x0^2 d12 − beta x0 x1 d14 + beta x0 x2 d24
    + (x0^2 − alpha x0 x3) d34
```

に同じ正規化を適用した。全alpha,betaについてS61の式を満たし、
`reports/S61_FAMILIES_20261030/almost_abelian_map.json` に61座標を保存する。
原点でのPfaffianは常に1なので、非退化性のための追加開条件はない。
一般alpha,betaで重み1成分が非零である。
v02のcompatible toric pencilの証明はMPS40の和への包含を与えるが、
今回の座標計算だけから特定のMPS個別IDを割り当てない。

## 検証の範囲

基底順序、Euler-gauge不変性、有限pushforwardの合成・逆、perspective正規化、
入力拒否・非Poisson負例は `code/verify_s61_adapter.py` で別に点検する。
134式の成功だけで新しい成分性や大域尽くしを主張しない。
一般族の写像、自然な一次部分への射影、残余fibreを別の記録として保持する。
