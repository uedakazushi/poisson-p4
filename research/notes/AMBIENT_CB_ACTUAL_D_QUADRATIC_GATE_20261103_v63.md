# Q_b側actual D：二normal blockには二次consequenceがない

2026-11-03。相対CB1∩CA=CI4の二normalを対象とし、CAの残る七方向の二次式を零と主張しない。

保存済みの全150×3 actual CB1−CA quadratic差を、**whole CB tangent38**へ射影した。
whole tangentは実二branch tangent37と真二次差Lで張られることが、先行の全ring/全38rankにより分かっている。
この150×38行列を物理的に再構成し、38次単元blockと逆行列で全150×3を除いた。
結果はQQ(b)上で恒等的に零である。
実kernel etaの二乗の係数rowを追加するとrank1となるので、このeta²はこのblockのactual D二次consequenceではない。

この射影がactual ideal全体に対する必要判定であることを確認する。
CB ideal元Fの一次項ellはwhole tangent38を殺す。
CB1とCAの二曲線は同じ一次方向を持つため、二次差の代入でFの共通Hessian項が消え、ell(raw quadratic差)だけが残る。
これは全38射影零により零である。
一次項が零のCB ideal元についても、その二次項はCB1の接空間上で零なので、CB1とCAの共通二normal上で零である。
従って単にいくつかのgeneratorを試したのではなく、actual Dの全二次idealをこの二normal planeへ制限したblockが零となる。
実GL×bによるintrinsic計算なので、smooth20 baseの一方向だけを固定した推測でもない。

先行diagramはZ4=R[u,v]/(u²,v²)→C2=R[v]/v²でuを殺す形である。
もし追加のambient equationがactual Dにu²相当の二次式を与えれば、nil moduleをC2へ縮める可能性があった。
この二次gateはそれを除外する。
しかし**cubic u²v**によるsocle killはまだ除外しない。
また七つの別方向を含む二次式f*z_iはこのblock外であり、whole Dの二次idealが零とは主張しない。

full D ideal、nil moduleの長さ3か4か、embedded primary componentは依然別の実cubic gateを要する。
非flat basechangeでfiber productのexactnessを仮定しない。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CB_actual_D_quadratic_gate_20261103_v63.py --output <new-output-directory>

completed98秒physical workerの全150raw二次差を再利用し、whole CBの全sourceと38次quotientを新しく計算する。
新しいGBやhigher-jet runなし、120秒・2.2GB。
