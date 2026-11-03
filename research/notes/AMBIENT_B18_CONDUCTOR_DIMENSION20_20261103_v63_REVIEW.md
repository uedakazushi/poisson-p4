# CA/CB原点conductor dim20・重複度上界4の独立検証

2026-11-03。対象は `AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63` のみ。
CAだけの後続multiplicity下界・full CM theoremはここでは承認しない。
新しい高次jetやGBは実行せず、既受理の全actual odd jetの入力hashを照合し、
今回の新しい接空間・IFT・初期idealへの論理橋を別に検証する。

## intrinsic conductor上の群作用

`A=R directsum I*t`、`B=R[t]/(t²−nu)` からconductorは `I directsum I*t`。
従って `C=A/conductor=R/I` である。個々のtarget成分はG不変、normalizationとconductorは
intrinsicなので、C上にG作用がある。選んだRがG不変であることは不要。

独立な七actual odd functionを消去して `C=S/L`、`S=C[[g19,W,U,V]]`、`L subset m²`
を得る。C上のLie derivationは各座標の像をSへ任意にliftしてSのderivationとでき、Lを保つ。
十九のconstant tangentを `partial_gi` とする。任意のh∈Lの二次初期式qについて、
`D_i h∈L subset m²` の一次部分は `partial_gi q` なので零。
従って全二次初期式はgに依存しない。既受理の全WUV restrictionから、その空間は正確に
`Q=span(WU,WV)` である。

hの全二次項を零にした後、三次初期式h3に対して `partial_gi h3∈(Q)`。
Qはgに依存しないので、g次数別にEuler恒等式を使えば、h3の正のg次数部分は(Q)に入る。
従ってg依存を二次式の倍数で除ける。formal product theoremや任意base changeとsource imageの
交換はこの初期式の議論に必要ない。

## 全physical接方向とactual七制約

元P0の二次monomialと外積の指数を直接微分し、全25個のphysical GL方向を作り直す。
producerのLie action helperはimportしない。独立のFraction Gaussian消去でG19、
G19+WUVのrank22、七補方向を加えた29×29 unit、G19+W/Vのrank20とunion rank21を照合する。

保存された全150 actual odd一次・二次・三次係数からh63の定数線形結合を作り直し、
全R-firstjetと全WUV二次項が零であることを確認する。上のderivation補題により、
actual七制約を解いたh63は全22変数で三次以上である。
exact monic U曲線で全actual odd値は三次以上、h63の三次係数は
CAで12/125、CBで−4/243。七制約を零へ解くIFT補正の先頭29成分を独立に再計算し、
七つすべての三次値を打ち消す行列恒等式を確認する。
補正はO(U³)、h63のR-firstjetは零なので、その値の変化はO(U⁴)。
またh63は最小表示で全三次以上だから、曲線の二次accelerationは三次評価に寄与しない。
従って真のcubic initialのU³係数が非零となる。

## 真の初期idealと比較

W/Vの両exact monic曲線を全150 critical式へ直接代入し、恒等的零を確認する。
Wのphysical曲線も元P0+W*x²∂z wedge ∂wへ照合する。
従ってg依存を除いたcubicのW³とV³係数は零で、WU/WVの倍数を除くと

    f3=U*(aU²+bUV+cV²), a!=0

を得る。b,cは未知のままである。
`(WU,WV,f3)` のHilbert–Burch matrixの全maximal minorを記号的に確認する。
heightは正確に2：height1 primeならWU,WVからWを含む必要があるが、Wはf3を割らず、
一方(U,V)はheight2の包含primeである。従ってHilbert–Burchでperfect/CM、
三変数商のHilbert seriesは `(1+T)^2/(1−T)`、重複度4。
g19延長からactual gr(C)へ全射があるのでdim C≤20となる。

G19×Wという実際のconductor内のsmooth familyのrank20 minorは非零である。
formal inverse theoremからそのsmooth20次元像を得てdim C≥20。
**この次元一致の後に**同次商のmultiplicity比較を使い、actual e(C)≤4と結論する。
G19×Vもsmooth20次元で、二接空間のunion rank21からW/Vの二つのminimal reduced support germは
異なる。追加のminimal branchやembedded構造を排除していない。

off-diagonal normalization pairはconductor上の有限fiber productの開部分に含まれるので
dimension≤20となる。この開部分自体がfiniteであることを仮定しない。
cone conductorはdim20/edim22、projectiveではdim19。
full I、actual CM/depth、multiplicity=4、CA/CB intersection Dは対象外である。
二つの異なる旗だけからgeneric conductor length≥2を推論することも行わない。

独立コードと全入力hash・証拠は
`code/ambient_B18_conductor_dimension20_review_20261103_v63.py` および
`reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63_REVIEW/` に保存する。
