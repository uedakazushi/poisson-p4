# rank一long残九源：独立同一点ambient監査

2026-11-08、PRIMARY E_multiple_hyperplane。**独立算術 PASS、B2/B4/B5/B6の全25次元成分性を受理。M1/M11/B1/B3/V1は未解決のまま。** 新規性は保存された既知成分との相違であり、外部文献全体の新規性調査とは別である。

## 1. 入力と独立性

producerの初回九点partial reportと、同じB6／V1だけのcontinuationを、complete-point index `1213b9368182a4c5cabeb41e845126b0bf7fa4ea56e0b9e096aa72d03d9409e3` のreport hash／case indexに従って結合した。元partialのoverall FALSE、B6の整数化guard停止と未完recordは変更しない。初回partialにはtop-levelのmonomial順序が無いので、明示されたindexの全150／350／126順序を独立標準順序と照合した。

B6だけ元tensor scale十二を二十四へ変えた同じprojective点であり、元の全150係数の定数二倍も直接確認した。元七点をproducerで再計算したとは記録しない。独立側は許可された同じ九点を一度ずつ再構成した。

独立callerは標準ライブラリのFractionと既存の独立sparse polynomial engineのみを使う。producer、SymPy、CASはimportしない。新しい点、minor探索、M1、profile又はjetは追加しない。実行前に29入力hashと180秒／3GiB契約のgateを確認した。

## 2. 元sourceからの全係数再構成

座標は(h,x,y,z,w)、r=y−x²/2、U=u∂z+v∂w、d=div U。保存150-vectorを点の定義として使わず、元g,U,H,Fから

```
Pxy=−F−(1+d)r, Pxz=u, Pxw=v,
Pyz=xu+H_w, Pyw=xv−H_z, Pzw=0
```

を構成した。全九点でg=0、(U+1)H=0、F=1+2z+3w+5z²+7zw+11w²を確認し、degree二斉次化と指定scaleにより全150係数を戻した。全strict10／projective350 Jacobi、全五Pfaffian、signed Euler Pfaffianの全126係数が一致した。実Dは `−scale²*homogenize_degree5(H)` で非零である。

元の全150方向を使った350×150 Jacobi微分D、全25 Euler方向、全25 GL方向とscaleを加えた150×51行列Mを独立に構成した。全entryの一致、既存nonzero minorの正確な有理数det、complete kernelの全元行消去とfree-coordinate identity、全D*M=0を確認した。新minorの探索はしていない。

| source | rank D | rank Euler | rank M | 射影接次元 | この点の軌道次元 |
| --- | ---: | ---: | ---: | ---: | ---: |
| M1 | 97 | 24 | 48 | 28 | 23 |
| M11 | 98 | 24 | 47 | 27 | 22 |
| B1 | 98 | 24 | 48 | 27 | 23 |
| B2 | 100 | 24 | 48 | 25 | 23 |
| B3 | 99 | 24 | 47 | 26 | 22 |
| B4 | 100 | 24 | 47 | 25 | 22 |
| B5 | 100 | 24 | 48 | 25 | 23 |
| B6 | 100 | 24 | 47 | 25 | 22 |
| V1 | 99 | 24 | 46 | 26 | 21 |

ONE実行は5.597267秒、maximum RSS 23356KiB、全九点PASS。report SHA256は `4b7c72ba014a53d85546b73626c357fc3388049f4d58cab64d561054009a1502`。

## 3. 四成分の認証と全境界

全七B1–B6/V1の純入口と別作者reviewは、各全sourceのinjectivity、正確な12→22→25次元、唯一の重複hyperplane、全frame／center／unit・zero境界を保持している。M1/M11の別作者純入口も同様に全25次元像を与える。これらの下界を数値点から推測していない。

特にB2/B4/B5/B6では、既約閉像Zの射影次元25が既証明である。Euler24方向を戻した元150係数空間のconeの局所次元は50以上であり、独立D100により全接次元は50。従ってこの点はregularで唯一の全成分を通る。Zはその全成分と同次元の閉既約部分なので等しい。全sourceとその全係数境界は同じ閉成分に含む。全点のsmoothness又は全点での局所唯一性は主張しない。

純入口の独立prior／相互比較を再利用する。B2/B5の残cubicの絶対既約性とmarked Hに対するsingular planeの位置、B4/B6のunique simple-L／rank三Qの制限rankは全族の不変量である。全prior129、非定数g新25、定数g二共鳴とscalar-cubic新三25への相違は既純reviewの全source証明に依存する。単に保存点のDが異なることや軌道次元だけで相違を結論しない。

M1/M11/B1/B3/V1については、正確な25次元閉像を通る点の接次元が28/27/27/26/26だったことだけを記録する。これらが成分、既知親の境界又は未知親に入るかは、本算術から決めない。

## 4. 一般軌道の純補足

`RANK_ONE_GZERO_VERTICAL_POISSON_FIELDS.md` を元全六括弧から別に確認した。任意affine ℓ(z,w)についてYℓ=ℓ∂yのLie微分の唯一残る成分は `−(U+1+d)ℓ`。従って全F／Hを保って同式の核がPoisson場となる。modular X=∂x+x∂yと独立な場は、M11/B4では∂y、B3ではz∂y、B6ではw∂y、V1では∂yとw∂yである。qの中心・frameを戻すときはこれらのactual pushを使う。

この全源stabilizer下界と同じ保存点のMの下界を合わせると、全九既約像のgeneric軌道は表の点軌道と一致する。四成分のgeneric軌道はB2/B5が23、B4/B6が22で、Rosenlicht商次元は二／三。これは全源のPoisson場を先に証明した結果であり、点Mだけからの一般化ではない。

元finite routing十七sourceでは、既閉八源と今回の四源を合わせ十二源が閉じ、未解決は五源となる。このsource数は成分総数と区別し、catalogue／global台帳は本reviewでは変更しない。
