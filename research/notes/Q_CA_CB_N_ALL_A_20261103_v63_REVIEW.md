# 全Qa source Smith・非零Qa完成環・actual交叉transportの独立監査

2026年11月3日。監査対象は manuscript_sync の三つの凍結packet
`Q_CA_CB_N_SOURCE_SMITH_20261103_v63`、
`Q_CA_CB_N_ALL_NONZERO_COMPLETED_MODELS_20261103_v63`、
`Q_CA_CB_ACTUAL_INTERSECTION_TRANSLATION_20261103_v63` である。
原source・report・既存独立receiptは変更しない。

## 1. 独立再構成と係数環

監査source `code/q_CA_CB_N_all_a_review_20261103_v63.py` は、対象producerをimportしない。
物理的raw CA/CBのtensor、divergenceによるcanonical射影、150係数の配列、
全25個のGL微分を別実装で再構成する。
三producerを新しいignored work先へ再実行し、runtime以外の全report項目と
Smith証明書の全byteも原reportと照合する。

物理座標 `(h,x,y,z,w)=(old0,old1,old4,old2,old3)` のpermutationを、
bivectorの符号と係数monomialの双方に適用する。
元の150係数で非零なのは
`90:a,105:-2/9,121:4/3,137:2/3,140:1` と厳密に一致する。
全CA20/CB20/N19 source方向とGL25方向の三行列を全項照合し、
各QQ[a] weight blockについて `U M V=D`、全非対角項の零、
左右determinantが非零有理定数であることを直接検査する。
CAとNに非定数Smith因子はなく、CBの唯一の非定数因子は `a/9`。
従ってsource微分像のaffine rankはCA28全a、CB28非零a/27零a、N27全aである。
これはtarget tangent dimensionの計算ではない。

## 2. 全非零物理Qaと同じS61曲線

物理的旧座標からの有限frameを

    Y0=x0, Y1=x1/a, Y2=x4/a², Y3=x2/a, Y4=−9x3/(2a)

とし、bivectorをさらに `a²` 倍する。
この変換のinverseに必要なparameterはaだけである。
全150canonical係数を、immutable S61再構成の曲線

    u25=−4/3, u36=3, u38=9a/2, 他0

と一致させる。原134式も全QQ[a]で零になる。
物理的a=0にはこのframeを使わず、S61曲線のa=0延長点と混同しない。

accepted closed proper flag-sourceを用い、上の全S61曲線を
固定β `(0,0,4/3,3/2,0,0)`、`r=0,k0=9a/2,他k=0` に同定する。
有限flag作用の後の全closed-source constraintsを計算し直し、CA/CB/Nの三つとも
QQ[r,z,w,a]のidealが厳密に `(r,z,w²)` になることを検査する。
ここでaは多項式変数であり、QQ(a)へのgeneric化はしていない。
定数βのintrinsic contraction kernelは `span(e_y,e_w)`。
chart外の唯一の候補e_wはa-independentなbase constraintで除外される。
従って全projective旗fiberが一点支持・scheme長2であり、未検査のa例外はない。

proper smooth birational sourceはこの零次元fiberの近傍で有限normalizationとなる。
有限normalizationのfiber長2から非正規性が従う。
単なるsource微分rankの減少を非正規性の証明に使っていない。

## 3. 全像環とfree k0

既存の独立監査済み全像定理は、七つのβ群因子を分離したsliceで
`A=R+I*p+p²*B`、`B=R[[p]]` を与える。
Rとp²,p³は実際のtarget座標で、p²BもAに入る。
従って全45個の一次odd係数が全image ringを決めるのであり、jetからの推測ではない。

監査では全61source mapからq0,q1,M0とactual selected-K block Lを作り直し、
全45項の `q1−M0 L^-1(q1)_selected` を独立に照合する。
その定数係数rowspaceは正確に

    N:     k1,k3,k4,k6,k7,k8,
    CA/CB: r,k1,k3,k4,k6,k7,k8

であり、sを加えたIのcodimensionは7,8,8である。
k0,k2,k5の列は全部零、pure-k0軸で全odd項は恒等的に零となる。
固定βをB=C=1へ戻すSp4 torusは
`B'=B*mu/lambda², C'=C/(lambda*mu²), k0'=k0/(lambda³*mu)`。
有限potentialを直接変換してこの式と他fiber係数が生じないことも確認する。
mu^5=27/16、lambda=(3/2)/mu²はaと無関係な非零複素定数として選べる。

従って全非零物理Qaで個々の被約成分の完成環は

    CA,CB: C[[u19,s1,...,s8,p²,p³,s1*p,...,s8*p]],
    N:     C[[u19,s1,...,s7,p²,p³,s1*p,...,s7*p]]

である。次元/embedding dimension/depth/multiplicityは28/37/20/2、
28/37/20/2、27/35/20/2となる。a=2/9も含み、全Poisson局所環の共鳴とは区別する。

## 4. 同じactual targetでの交叉transport

個々の完成環の同型だけからCAとCBの実交叉の同型は従わない。
ここでは両全source mapに同時に可換するtargetの有限automorphismを使う。
source座標順 `(r,k0,k1,k2,k3,k4,k5,k6,k7,k8)` において

    (k0,k2,k5) -> (k0,k2,k5)+t*(1,−p,p²/3)

を両族の全61多項式に代入すると、厳密にu38だけがt増え、他60項は不変である。
actual selected-K blockのdetは両族とも恒等的に−72。
同じ有限置換のtを−tにしてinverseも直接照合する。
従ってこのtarget translationは両全image kernelを同時に保ち、IA+IBも保つ。

さらに独立監査済み49変数actual presentationのCA109式とCB109式を全部読むと、
全218式が文字通り `v4=u38` を含まない。
nilpotent F、保存されたpFのtarget expression G、colon primeの全45生成元も
v4を含まない。ゆえに古いexact ideal/colon/nilpotence receiptを保ったまま、
全actual intersection環とそのnilmoduleの作用をfree-v4方向へ運べる。
normalization同士のfiber productをactual交叉の代わりに使っていない。

固定βtorusもaに依存せず、全Qa曲線はv4だけの直線に移る。
従って全非零物理Qaのactual DはQ1と同じ完成環を持つ。
S61ではdim18/depth11、非零主nilradical(F)はsquare-zero、nilmoduleは
滑らかな11次元primeの上のrank1、associated primeは正確に二つである。
九つの値・perspective因子を戻すとdim27/depth20、nilmodule支持は20次元となる。
環の直和分裂は主張しない。

## 5. a=0とscope

物理的Qa0には全非零Qaの完成型・fiber長・actual Dの型を延長しない。
個々の被約成分についてのみ、excellent varietyの非正規locusの閉性から
Qa0も非正規である。既存の全射影R1定理と合わせるとCMでもない。
これはQa0でclosed-point depthが2未満という主張ではない。
CM局所環ならその全localizationsはS2であり、R1と合わせてnormalとなるという議論である。

全Poisson空間の枝の尽くし上げ、a=0の全完成型、generic-rank2の全global交叉境界は対象外。

## 再実行

`code/q_CA_CB_N_all_a_review_20261103_v63.py` に三つの新規producer replayディレクトリと
`--output reports/Q_CA_CB_N_ALL_A_20261103_v63_REVIEW.json` を渡す。
最終JSONが全検査の成否と全入力SHA256を保存する。原frozen packetに変更はない。
