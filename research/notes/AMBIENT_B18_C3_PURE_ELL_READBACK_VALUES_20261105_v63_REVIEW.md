# 既存 pure-ell164 と physical12 の独立 VALUE 検算

2026-11-05。PRIMARY report
`AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/report.json`
SHA256 `341251643c4936ca11d818ee69f7c8a336898cbea86b3bebd4b6dc8ad825d208`
と、準備・SOURCE peer `58abe686…` を固定して検算した。
実行済み PRIMARY は strict temporary-height guard を持つ V2
`29308014c69da649c66596454c95a749c14d8724f129c4ac97f9cfb939da1ca4` である。
元V1と全source bytesは変更されていない。

独立 REVIEW report:
`reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json`
SHA256 `2c2c157f217212c14c7cc46a94f252701d9b6909e1c3f13fdc5b470351c9ecd1`。

同じ provided105 の315係数を入力に、元18584演算を別の jet2 実装で評価した。
author helper/producer/native code は importしていない。
畳込みは係数・summandとも逆順、unit inverseは直接の二次公式ではなく

\[
x=1+u\quad\Longrightarrow\quad x^{-1}=1-u+u^2\pmod{\ell^3}
\]

として構成した。全中間値のstream hashを保存した上で、以下を確認した。

- 全164 root の三係数が保存された PRIMARY 値と一致する。
- 元105式が実際の F9 target (z_A=r_A=0,s_A=\ell/5) を満たす。
- 全164 root の定数・ℓ一次係数が受理済み C2 cache と一致する。
- 六つの inverse すべての逆順の積が1となり、二つのνも元Q20成分と一致する。
- 同じzero-based physical150の12行について、aA−aB、bA、bBの全三係数が一致する。

結果は pure-ℓ の二次係数に関して

\[
\delta_{\ell^2}=0,\qquad b_{A,\ell^2}=0\quad(\text{全12行}),
\]

\[
b_{B,\ell^2}(18)=-1/15,\qquad
b_{B,\ell^2}(60)=-29/120
\]

であり、bBの他の10行は零である。値は検算前に零と置かず、元DAGから求めた。

Fraction の各和・積は、分子・分母のbit数から演算前に上界を確認した。
保存DAG係数の観測最大18 bit、一時QQ結果の観測最大18 bit、
演算前上界の観測最大30 bitで、2048上限内に収まった。
非自明な guarded QQ演算2013回、実行時間約0.253秒、
入力・metadata・値・scratch込みの静的上界128327で130000未満、
RAW入力＋出力477174 bytesで2 MiB未満である。

これは既存キャッシュの固定 readback の検算であり、新しい陰関数、gauge、source jet、
C3のrank、三次 closure、torsion指数の判定ではない。
