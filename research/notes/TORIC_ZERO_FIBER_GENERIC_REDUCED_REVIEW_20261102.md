# T01 零ファイバーの一般被約性：独立レビュー

2026年11月2日。[主証明](TORIC_ZERO_FIBER_GENERIC_REDUCED_20261102.md)を受理する。
指定 AB 点で T01 は滑らか16次元、実際のスキーム論的零ファイバーは滑らか10次元である。
従って AB の一般点での零ファイバーの長さは1。
全 AB 上の被約性や、全 S61 の局所構造の主張ではない。

## 1. 元134式からの別 Taylor 計算

独立コードは主 producer を import せず、legacy の134式から再構成する。
各単項式の変数重複度を Counter で集約し、二項係数を使って Taylor 係数を求める。
主 producer の occurrence ごとの微分ループとは別の実装である。
点は既存の対数退化証明書の S61 座標から再評価し、保存点と比較する。

元 Jacobian の全階数43、一次列の階数10、二次列の階数33、および43次 minor を確認する。
保存18方向は別に rank18 と J U=0 を検証する。
陰関数消去は row と column の pivot 候補をともに逆順にして取り直し、
二次残余の全134行を計算する。得られた行空間は主の七行と厳密に一致した。

共通 Hessian kernel12次元、活性六変数、七行の可逆変換も確認し、

    J2=(u,v)²+(ua,ub,vc,vd)

を得る。変数変換の determinant は8。
この単項式 ideal の商を別に次数0,…,8で数え、
k[a,b,c,d]⊕u k[c,d]⊕v k[a,b] の寸法公式と照合する。
12自由変数を加えた下側初期商は dim16、multiplicity1、radical=(u,v) である。
これは真の初期 ideal に含まれる七式であり、全初期 ideal との等号は不要。

## 2. multiplicity1 の適用範囲

既知の被約既約16次元 T01 は指定点を通る。
従って下側初期商から gr O_T01 への全射があり、両者の次元は16。
Hilbert 関数の比較から e(T01)≤1、非零局所環なので e=1 である。

T01 の局所環は複素有限型整域の局所化で excellent。
completion は被約かつ等次元なので unmixed である。
この範囲では multiplicity1 による正則性判定を適用できる。
主ノートの有限 Noether normalization による説明も妥当である。
一般線形 parameter は maximal ideal の reduction となり、その正則冪級数環上の
有限 torsion-free rank は multiplicity1。有限整拡大は rank1 なので分数体に埋まり、
整閉性から正則冪級数環自身に一致する。

参照した
[Huneke–Swanson の著者公開本文](https://www.math.purdue.edu/~iswanso/book/SwansonHuneke.pdf)
Example 11.1.11 と Exercise 11.8 はこの multiplicity 判定を記載し、Nagata に帰属させる。
本研究で使うのは上記の excellent・被約・等次元 completion の場合である。
周囲 S61 の embedded 部分を含む環に同じ判定を適用してはいない。

## 3. β の零切断

正則16次元となった T01 の接錐は線形 prime である。
J2 を含むこの prime と radical(J2) はともに height2 なので等しい。
従って T01 の接空間は元18方向の中の u=v=0 である。

独立計算でこの16次元空間の β 微分階数は6、保存 minor は −3/4。
六つの独立な β 座標で切った局所環は正則10次元となる。
その中に既知の10次元 AB が閉包含し、全 β は AB 上で消える。
正則局所整域の非零 ideal では次元が落ちるため、この閉包含の ideal は零。
従って残りの十個の β も消え、全16個で切った実際の零ファイバーは同じ正則局所環である。
ここでは full S61 の18次元接空間から直接 smoothness を結論していない。

既知の T01 の次元・被約既約性、指定点の対数退化による所属、
零ファイバーの被約台 AB という等号は既存証明への依存である。
元134 Taylor2・七式・微分 minor は本レビューで独立再構成した。

    .venv/bin/python code/toric_zero_fiber_generic_reduced_review_20261102.py --output reports/TORIC_ZERO_FIBER_GENERIC_REDUCED_REVIEW_20261102.json
