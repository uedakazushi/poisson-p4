# 全conductorと任意nil-cycleによる上界7：独立監査

2026-11-04。元の actual CA/CB 交叉のB18代表点完成環において
`Nil^7=0`、指数範囲 `[3,7]` を受理する。旧上界12の frozen 資料は変更しない。

## 同じsource上の11節点 trace gate

元41315節点 source と joint105 のprefixが一致すること、全105 bindingと
キャッシュの元変数順、単項式順 `(1,r,s,...)`、同じ physical R0 九行を確認する。
元のνA=20086、νB=40865は、正規化chainの同じ trace-free monicの根でもある。
既存105 jetを使い、依存11節点だけを独自の前向き一次演算で再計算すると、
両方とも `(constant,dr,ds)=(0,0,0)`。

独立監査済みの exact `J=IA+IB=(zA7,r²,rs,s²)` は `m_R0²⊂J` を与える。
従って両νはJに入る。これは関数が恒等的に0という推論ではない。
既知の二乗最大idealを含むidealへの所属なので、一次値で全次数の判定が終わる。

## すべてのnil-cycleのeven係数

actual free4環は同じ原点の
`C=R0[τA,τB]/(τA²−νA,τB²−νB)`。
受理済み centered proper/physical interface を固定し、native dとR9各関数が
原点最大idealに入ること、および全36 Qrad列のC要素が `m_C²` に入ることを
源構造から追う。basis倍もこのidealを保つ。

4 basisの積でeven成分が出るのは同じmask同士だけで、係数は順に
`1,νA,νB,νAνB`。よって `f,g∈m_C` なら
`(fg)0=f0g0+νA fA gA+νB fB gB+νAνB fAB gAB∈J`。
従ってQradのeven行全36係数はJに入る。

完全H1≅Nil(D0)により、任意のnil元はcycleから回復した
`b0+Σci βiτA` の形を持ち、`b0=Qeven*a∈J`。
β8は実 `(zA7,gA)` であり、旧gAを落とさない。
特定の候補cycleや未証明のnil生成元だけを調べた議論ではない。

## 全conductorによる基底変換核の消去

二つの正規化を行う順序を入れ替えても、各第2段階の任意基底変換核は平方零。
従って合成核Kは `K²⊂KA∩KB` を満たす。
`IiBi⊂Ai` の任意のcについて、乗法 `Bi→Ai` はAi-linearであり、
包含との合成はAi上のc倍である。従って任意の基底変換後の包含核はcで消える。
有限な共通R0上のtensorは既に完備なので完成でこの議論を失わない。

even conductorだけからは `J K²=0` が従う。
さらに **全conductor** は各actual odd `βiτA` も含むので、これらもKA、従ってK²を殺す。
任意nil元の上の表示により `Nil(D0) K²=0`。
Torを捨てた単射や、scalar conductorだけではこのodd項の結論は出ない。

旧 `K^4=0` からK⊂Nil、従って追加的にK³=0。
whole E4のnil三乗零より `Nil(D0)³⊂K`、ゆえに
`Nil(D0)^7=Nil(D0)*(Nil(D0)³)^2=0`。
E4のembedded N/Pでnil二乗が0とは仮定しない。

全正次数nil積の消失、完備化での閉性、actual G19およびscalar chartの
既受理gateを使うと全B18実germへ上界が移る。Qbのcoherent nil-square閉包による
下界3は保持される。正確な指数、三乗の零・非零、H0全体のannihilatorは未決。

再現：`.venv/bin/python code/ambient_B18_quadratic_conductor_nil_upper7_review_20261104_v63.py`
