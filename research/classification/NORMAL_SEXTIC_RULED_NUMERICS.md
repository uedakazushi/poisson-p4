# normal degree六のirregular ruled残枝：下降・fiber次数・nef数値還元

2026-11-07。主義務G。出発点は独立確認済み
`RANK_ZERO_NORMAL_SEXTIC_REDUCTION.md` とその separate review。
S⊂P4はnormal・integral・非退化degree六、rank零、実Dは既約被約quinticとする。
既にq(S)>0、sectional genus g≤3、全特異点rational、解消Xはirregular ruledのblowupと分かっている。
smooth Sは既受理のsmooth非plane定理で除外済みなので、本稿ではSがsingularな場合を絞る。

## 1. rulingの下降とfiber次数一の排除

ρ:X→Sをminimal resolution、p:X→Bをminimal ruled modelから引き戻したrulingとし、b=g(B)=q(S)>0。
ρのexceptional curveはrational singularityの解消の連結rational treeである。
正種数BへのP1からのmorphismはconstantなので、pはρの全fiber上でconstant。
(ρ,p)のreduced image Γ⊂S×BはS上proper quasi-finite、従ってfiniteでbirational。
SのnormalityからΓ→Sは同型となり、g_S:S→Bが下降する。

S→Bは各base DVR上torsion-freeなのでflatである。normal surfaceはCMなので各fiberもCMの純一次元。
generic fiberはP1で、H=O_S(1)のfiber次数をa≥1とすると、全fiberのHilbert polynomialはa*t+1。

a=1なら、ample Hとdegree一から各fiberのfundamental cycleは唯一のmultiplicity一のline。
CMとgeneric reduced性からembedded point又は零次元nilpotent構造は無く、fiberはschemeとしてP1。
全fiberがsmoothなのでflat morphism S→B自体がsmoothとなり、Sもsmooth。
従ってsingular残枝では厳密にa≥2である。

一般hyperplane Cは有限個の特異点を避け、滑らかなgenus g、degree六の曲線となる。
C→Bはdegree aの有限写像なので b≤g 及び `g−1≥a(b−1)`。
g≤3ではO_C(H)はnonspecialでh0=7−g。
linearly normalなSのhyperplane完全列から、dimension3−gのcokernelがH1(O_S)へ入り、b≥3−g。
従ってg=0,1は不可能、残るのは

`(g,b)=(2,1),(3,1),(3,2;a=2)`。

## 2. singular Sにsmooth double-point formulaを用いないnef恒等式

ε:X→X0をB上のminimal ruled surfaceへのblowdownとする。
total-transform exceptional basis E_iを取り、L=ρ*Hを

`L=ε*L0−Σ m_i E_i`

と書く。Lはnefで、各total E_iはeffectiveなのでm_i=L·E_i≥0。
各blowup段階で、その中心のtotal fiberから新exceptionalを引いたdivisorはeffective。
後続blowupへpullbackしてもeffectiveなので `a−m_i=L·(F−E_i)≥0`。
infinitely near中心もこの議論に含まれる。
**必要な境界は0≤m_i≤aであり、smooth ampleの場合のstrict m_i≤a−1ではない。**

交点式は

```
L0²=6+Σm_i²,
K_X·L=−(6+Σm_i²)/a+2a(b−1)+Σm_i=2g−8.
```

最後の等号はcanonical Weil classとCartier Hのprojection formulaによる。
従って全blowup数・全infinitely-near配置について

`6=a(8−2g)+2a²(b−1)+Σ m_i(a−m_i)`。                 (★)

この式にはSのChern数又はsmooth double-point formulaを使っていない。
全summandは非負なので、(g,b)=(2,1),a≥2は右辺≥8で矛盾。
(g,b,a)=(3,2,2)も右辺≥12で矛盾。

## 3. 正確に残る二つの数値型

残りは b=1,g=3のみで、(★)は `6=2a+Σm_i(a−m_i)`。
従ってa≤3であり、singular残枝a≥2から

| a | 残るmultiplicity |
| --- | --- |
| 3 | 全m_iは0又は3 |
| 2 | m_i=1が厳密に二個、他は0又は2 |

と全数値支持が尽くされる。m_i=0/aの段階を非正則なelementary transformで除いたとは扱わない。
有限の交点恒等式と整数支持確認は別のFraction証明書に保存する。
ここまでの還元自体は、後続のrelative-conic又はelementary-transform分類を仮定しない。
