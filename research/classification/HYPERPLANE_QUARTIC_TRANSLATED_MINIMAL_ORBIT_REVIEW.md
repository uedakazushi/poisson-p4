# 最小冪零・長block translationの三軌道と一般成分からの排除

2026-11-07。主義務G。対象は、単純Hと既約四次Fを持ち、Hのmodular cofactorが
Jordan(2,1,1)、translation top長2となる型10の全支持である。
`HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS.md` の正常形論を独立読解し、
`code/classification_translated_minimal_orbits_review.py` で主要恒等式と全三軌道rankを再構成した。
証明書は `CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_ORBITS_REVIEW_20261107.json`。
producer/CASをimportせず、既存の独立Fraction・疎多項式算術だけを用いた。
実行0.192秒、最大RSS14,720 KiB、メモリ上限1 GiB。

## 1. 正常形完全性の確認

Z=y∂x+∂yを固定した元の全90係数源に対し、別に25自由parameterの全六括弧を生成した。
25×25のaffine変数変換行列と保存された逆行列の積を両順で単位行列と確認し、
全90係数が既存完全源と一致した。元の次数二制約、divergence、cofactor不変性も確認した。
したがって補助座標を先に課した小切片を全支持にすり替えていない。

補助(r,s,z,w)座標の全JacobiがΨ=VΦ+Φの三偏微分に等しいことを再計算した。
δ=0なら最高四次Pfが零となり単純Hに反する。
δ非零でα=1/2、`δA₂=L/2+Vbase(L)` が従い、
Y=δr+L、Htilde=2δG−L²と置いた**分母なし**の完成平方恒等式を確認した。
その結果が元座標で `x_new=x+L(z,w)/δ` という可逆affine変換であることを独立読解した。
r=x−y²/2自体をPGL変換には使用していない。
potential定数を引く操作はtensorを変えない。

trace−1/2の二変数affine場について、次数三以下の全単項式を独立に列挙した。
weight−1が現れるのは、座標交換を除いて
`(-1,1/2):z`、`(-1/3,-1/6):z³`、`(-1/2,0):z²,z²w` の三型のみ。
重複固有値−1/4は次数三以下にweight−1を持たず除外される。
最後の型の零固有値translationも残し、既約性がb非零を強制するのでそのtranslationは零。
残る平行移動とz,wの非零scaleは全てaffine/PGLであり、z³型で必要なcube rootも
複素数体上の有限な非零代数選択にすぎない。
よって全対象点は、A/B/Cの**高々三つ**の実affine/PGL軌道へ入る。
三軌道の相互識別は今回の弱い結論の追加前提ではない。

## 2. 全GL5・垂直・scale行列

三代表を元五変数の斉次二次bivectorに戻し、保存された整数witnessと全係数一致を確認した。
各代表で全25個のGL5微分、全25個のEuler-vertical微分、scalar方向を別に生成した。
各行列の全有理数rankと非零minorを再構成して、三者とも

    vertical rank=24,
    vertical + scalar rank=25,
    vertical + scalar + full GL5 rank=47

を得た。従って実際の射影PGL軌道次元は `47-25=22` である。
接空間の大きさを軌道次元と取り違えていない。
A/Bの全空間成分性や成分帰属は、この軌道rankからは結論しない。
Cの既知MPS帰属は別の全重みslice・実半単純曲線の独立証明に依存し、本補題には不要である。

## 3. Gへ渡す最弱の結論

型10の対象全支持のPGL飽和は、高々三つの22次元軌道の和である。
従ってその閉包の射影次元は22以下。
既存の `DOMINANT_COMPONENT_GEOMETRY.md` による未知G成分の次元25,26,27という制限と
合わせると、それらの成分の一般元がこのcofactor型になることはない。

これは**22次元stratum自身のprimitive β像が非優越**という主張ではない。
22次元空間から16次元底への写像は優越し得る。
一般β上の全点の被覆、A/Bの例外支持への帰属、他のcofactor型の排除も未解決のまま残す。
新成分数の追加は行わない。
