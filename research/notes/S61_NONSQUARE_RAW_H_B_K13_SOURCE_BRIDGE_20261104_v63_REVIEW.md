# 有限 Ka13 source bridge の独立検証

対象は凍結済み MIXED_MINOR_W、K13_E_RECONSTRUCTION、K13_H4_INPUT、
K13_MIXED23_SCALARS の四つである。作者の worker や arithmetic helpers を import
せず、別実装の Python 整数係数積和・整数長除法・分母を払った交差積だけで検査する。
新しい solve、逆元探索、polynomial GCD、norm、点、有限体、GB は実行しない。

元の adjugate column への一つの共通 QQ*a^24 normalization を保持し、その16成分と
mixed norm matrix の全256係数の g13 reduction を独立に source へ結び戻す。全16の
W raw products と g13 divisibility、V15 の literal unit Bézout、および全16の
companion identities EV=tV を検証する。

H4 は元の Laurent74項から五つの係数を作り、a の閉じた逆元、全ての保存
pseudodivision、正規化を確認する。最後の係数が1であることも独立に調べる。
m2,m3 の元の pencil coefficients への substitution と零性、prepared B7×3 の
345 source terms からの全84係数（69非零・計789 QQ 項）を確認する。
既存入力・証明書は変更しない。

有限射影の結論は次である。T=QQ[a,a^-1,e]/(h1,k1,m1) の QQ length は、既に
独立検証された共通 monic16 presentation と Fitting/Smith の議論で13以下である。
元の adjugate relation と W=0 から k1(E)V=m1(E)V=0。companion identity
h1(E)=0、EV=tV と V15 が単元であることから h1(t)=k1(t)=m1(t)=0 が従う。
従って T は QQ[a]/g13（QQ dimension13）へ全射し、length bound により同型である。
g13 の既約性・reducedness は不要。monic H4 を加えた containing source cover の
QQ length は52である。

これは actual Bbad の存在・空性、coker B の module length、全 cubic57 や
quartic239 の証明ではない。D(G) での actual source membership と G0 の扱いは、
それぞれの既存証明を別途保持する。実際の検査結果・入力 hash・時間は report と
independent-proof に記録する。
