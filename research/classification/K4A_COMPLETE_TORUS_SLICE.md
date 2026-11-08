# K4a の全 equivariant 切片は次数五で完全：最弱の成分性/実親二択

2026-11-08、PRIMARY E。これは既Nの高次jetを際限なく延ばす経路ではない。
全torus表現が exact formal graph の次数を事前に五へ制限するため、
一つの有限 polynomial family と全350式の全次数十を確認すれば全切片の被約支持を判定できる。
pure経路は root と generic 作者が別に監査済み。計算の独立再生は別に要する。

## 入力と true 102 切片

P は凍結済みK4aの同じ18倍原始整数点。全Dは350×150、rank101。
実Euler/GL5/scalar tangent像Gはrank48、Eulerはrank24、kerDは49次元。
保存kernel列23の N は G 外、全Gを殺すfunctionalで値8/87を持つ。
新しい point、minor、rank、normal complement を選ばない。

既存D101証明書の独立係数列を U=(E_j)_{j∈C}、独立101方程式行を R とする。
D_RU は非零QQ determinantを持つ同じ保存blockである。
各 E_j は一つのraw二次係数なのでtorus eigenvector。

```
S=P+t*N+sum_j x_j*E_j
```

は102次元actual線形切片。D Uが可逆、N∉G、G⊂kerDから
`A150=G ⊕ span(N,U)` が従う。新150minorを探す必要はない。

## Equivariance が full IFT graph を有限化する

normalized作用は `s^-1 diag(s²,s,s,1,1)_*`、Pを厳密に固定する。
全150raw係数のweightsは[-5,2]、Nのcharacterは-1。
各元350式のweightは quadruple出力重みの和−quartic monomial重み−2。
保存Dの非零entryは同じinput/outputweight間にしか無いことを全supportで確認する。
従ってSと選択101方程式もT安定である。

選択101方程式の∂/∂x=D_RUが可逆なので、formal IFTは一意なgraph x=f(t)を与える。
f(0)=0、DN=0よりf'(0)=0。
一意性からgraphはT-equivariantであり、weight w_jのcoordinateには

```
f_j(t)=constant*t^(-w_j)
```

しか現れない。w_j≥0では零、w_j=-1でも一次係数零。
従って全formal graphは exact polynomial

```
P_t=P+tN+t²X2+t³X3+t⁴X4+t⁵X5
```

で尽くされる。XjはUのweight -j部分だけを持つ。
次数二から五の固定block再帰は、この既に有限と証明された四係数を決める手段であり、
未証明のpower-series remainderを残すjet実験ではない。

## 全350の degree≤10 判定

元Jacobiはraw係数に二次なので、P_t代入後の全350式はt次数高々十。
選択101式も含め、全degree0,...,10の係数を厳密に確認する。
選択101式が全次数で零であることはexact graph証明の必須確認で、次数五で止めない。

全350が零ならliteral actual familyが得られ、formal切片はsmooth A1。
Gの48個の実作用parameterを取った action×S 写像は原点で微分が同型、従ってétaleな局所積。
閉軌道を仮定する Luna theorem は使わない。Poisson方程式は作用に不変なので、
全raw局所schemeはsmooth49、Euler24とscale1を除いてsmooth unique projective24である。

逆にfull350のどれかに非零係数があれば、各式はtorus weightから単一monomial c*t^k。
その一式だけで全formal切片の被約支持はt=0。上の実action局所積からraw被約local dimension48、
projective dimension23となり、既actual23orbitの閉包が全ambient被約成分である。
schemeのsmoothness、reducedness、正確なnil指数を主張しない。

どちらの結論も全同次切片の確認であり、既P+eN pencilの非零自己Jacobiだけからは得られない。
追加point・normal・次数増加は行わない。新成分の計数は独立完全再生と全prior126相違監査後の統合担当のみ。
