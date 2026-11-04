# R434とactual Saux：固定商剰余の独立検査とgeneric chart

2026-11-04。同じF32003[a,e,d]の保存商3212項・剰余6729項を読み、
新しい除算・factor・点・素数なしに `qR+rem=Saux` を全係数で検算する。
DEGREVLEXのR leadingは `(2,15,3)`、係数27306。
全6729剰余単項式がこのleading monomialで割れず、剰余は非零である。
従ってRbarはSauxbarを割らない。

## 実sourceから標数0へ

先行する独立記号Fp checkerは、768個の実source演算と公開7022係数の還元が一致することを
確認済みである。今回さらに同じsource-bound入力の全QQ分母のp-unit性と、
保守的a冪下界を追う。最終下界は非負なので、実際のcleared source式Sauxは
Z_(32003)[a,e,d]の多項式である。公開QQ展開の厳密全係数を未独立監査のまま信じて
この多項式性を推測する必要はない。

Rの原始整数模型はp-primitive。もしQQでRがactual source Sauxを割れば、
Gaussの補題により商もp整となり、還元後にも割るはずである。
非零normal remainderがこれと矛盾する。
別の独立F32003既約性gateからQQでのR既約性があるため、
`gcd_QQ(R,actualSaux)=1`。

RとSauxはQQ多項式環で正則列となる。平坦なCへの体拡大は、Rで割った商上の
Saux倍の単射性を保つので、Cでも共通divisorを持たない。
R自体が絶対既約という仮定は要らない。
これは多変数idealがunit、または零点集合が交わらないという意味ではない。

## 他の分母と実P5 coverage

源に束縛されたa次数は、Rが11、pが2、Delta407が8、h22が4、Bが1。
後四つは非零でRより低次数なので、QQ既約Rはどれも割らない。
受理済み恒等式 `Delta5=-(6/5)B p² Delta407` より
`gcd_QQ(R,Delta5*h22)=1`。全FULLU named factorsもa次数が1以下である。
これらのcoprimalityも平坦な体拡大で保持される。

H5だけを最初に消去するD(Delta5)上では、u4のheavy行はpair-only。
実補正pure行は `rho_u/(delta2*rhoLCM)*Bhat` であり、P4にその行とtarget15を
付加したP5 determinantは `rho_u/(delta2*rhoLCM)*Saux`。
この式はD4を逆にしないborder determinant恒等式である。
従って `FULLU∩D(Delta5*Saux)` では実P5＋H5の10pivot chartが使える。

RとDelta5*SauxがCで共通divisorを持たないので、このchartは幾何的なRの各成分で
稠密である。R側で残る交わりは、非空ならcodimension少なくとも2。
残行列の階数、global cubic57、chart補集合の全除去、unitideal/disjointnessは主張しない。

再現：`.venv/bin/python code/s61_nonsquare_prePOWER_R434_Saux_divrem_review_20261104_v63.py`
