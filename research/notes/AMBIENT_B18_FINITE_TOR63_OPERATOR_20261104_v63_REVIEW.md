# Finite Tor63 literal operatorの独立監査

2026-11-04。対象は凍結された `AMBIENT_B18_FINITE_TOR63_OPERATOR_20261104_v63`。抽象Tor2、長さ63・最小27生成、自然な有限像と原点二次gateについては、既受理の `FINITE_TOR63_IMAGE_20261104_v63_REVIEW` を明示的に使用する。今回新たに確認するのは、その全次数regular-function operatorのsource bindingである。係数展開・非所属判定・高次jetは行わない。

## Unit105 cloneと90正則商

全209337 scalarノードの構文、原点値、unit逆元を独立lookup interpreterで確認する。元actual source全prefixと、九clone完成時checkpointがそのまま保持されることを比較する。元105式・105未知数のJacobianを独立疎ADで計算し、保存Jacobianと全係数一致、保存逆行列との両積、非零det=1202315964973056/390625を確認する。

九つのcloneは七prefix（長さ0,...,6）と二axisである。全105未知数の新しい名前・相互disjoint性、全source祖先のhomomorphic clone、全105式と九targetをliteralに照合する。target以外に元の未clone変数が残らないため、clone系のJacobianは同じ元unit105である。新しい求解を行わず、pointed unit解の一意性によりregular functionsとして定まる。

九函数それぞれの七prefix差は、対応するY_jを零にすると前後のtargetとunit系が一致する。従って63商は正則である。最後のprefix長7は同じtargetを持つ元source解なので、その一致も一意性による。各numerator/denominator・前後root・cached一次原点値を確認する。

残る27商は各plane函数のr²,rs,s²分解である。CA plane函数のorder≥3、全八CB plane函数のorder≥2を既受理cacheと照合する。axis値はそれぞれr²,s²に可除、mixed差は両axisで零であるため、正則座標環で(r)∩(s)=(rs)を使える。全numerator/denominator・cached二次原点値を確認する。これらは正則函数としての商であり、その分母の零上へ0/0を新しく代入しない。

## 実行列から全512逆恒等式

九coordinate telescopeと九plane分解から、元16生成元Gとcanonical10函数hに対するG=A*hをsourceへ結ぶ。全A成分は90個の新正則商とQQ定数だけである。selected10とredundant6は前preflightと同じで、quadratic3 matrix Uの原点detは625/768である。新しく逆にするregular商を含む函数はdet Uだけで、原点unitを確認する。

実DAGの90商を独立な自由変数へ置き換える。実det Uの根が通常の3×3 determinantに等しいこと、実Uiの九成分がadjugate/detに等しいことを確認する。実16×16変換Bと実逆行列Cをその回路から直接読み、Bの各成分に含まれる(det U)^-1が高々一つであることを確認する。det Uを掛けて分母を払い、自由多項式環で

`(det U * B)*C = det U * I16`、
`C*(det U * B) = det U * I16`

の全512 entryを独立に検算する。著者のuniversal checkerを再実行せず、実保存行列を起点とする。これによりB*Aの最初10行がI10、最後6行が零となり、canonical recoveryと六actual syzygyが従う。巨大なsource多変数係数を展開した証明ではない。

## 八cycleと二十七actual C0 map

六cycleの原点を前preflightと全係数比較する。二etaはcanonical三quadricの(-s,r,0)、(0,-s,r)を実Tで持ち上げたものとしてliteralに確認する。canonical10函数が(z7,r²,rs,s²)であることも束縛する。

joint105のCA blockはCBのslave blockから独立なpureCA57であり、同じF9 targetを持つ。unit解一意性により、CA plane f63は旧gA、七coordinate商はC0 free9で用いる正則係数と一致する。この全次数の同定を使用し、原点の一致だけからbasisを同一視しない。C0へのf_i tauA vector、八つの実even差、同じphysical row順を確認する。

十五q∧qと十二eta∧qについて、全h8×8と九成分imageの全根をlookup演算で別に再構成する。CB outputとの同一視はg_j tauB=(aAj−aBj)+f_j tauAを使い、乗法はCAのtauA²=nuAを用いる。二族traceを直接同一視しない。old r/gを再評価する操作は含まれない。

これらは既受理Tor2の自然像であるため有限かつm_R0²-annihilatedである。しかし像が零・非零であること、640関係への所属、socle、D0または全Dのdepthは今回判定しない。全input/evidence hash・時間・メモリは専用REVIEWへ保存する。
