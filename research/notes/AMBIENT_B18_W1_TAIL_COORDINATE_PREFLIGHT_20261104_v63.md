# Actual w=1 tailからS61への有限座標mapと、17方向の単位切断

2026-11-04。これは全218式とのactual germ bindingを構成するための有限入力preflight。
旧局所型を同一視してJinfのsyzygyを得た、という主張ではない。

## 1. 実際の150出力をw=1へ戻す

既凍結のCA sourceは各physical rowのw次数0..3係数を、各自のtrace basis
`(1,tau_A), tau_A²=nu_A`で保存している。その四係数を足して全150個の
`H_A=R0+m_R0*tau_A`値を持つw=1出力を作った。
原点値は既存のE(1,0)そのもので、tau係数の原点は全て0。
これは新しいbeta点ではなく、保存済みの同じphysical W1点である。

旧W1認識の固定GL5を適用した。固定GLがvalueをJにするのは原点だけなので、
近傍ではさらに実際の6個のvalue coefficient(a,b,c,d,e,f)に対するDarboux frameを入れた。

```
delta = a*f-b*e+c*d,
F4 = [[1,0,0,0],
      [0,1/a,0,0],
      [d/a,-b/a,1,0],
      [e/delta,-c/delta,0,a/delta]].
```

原点でa=delta=1。これは旧normalizerと同じ普遍的な有理恒等式
`F4*C*transpose(F4)=J4`を満たすframeであり、明示逆行列も保存した。
二次trace代数における逆元はnormのunit inverseで計算し、tau自体を逆にしない。
全係数は同じactual H_Aに留まる。

その後、既存Euler quotientと有限perspective変換を適用した。
全150正規化係数、全61座標、primitive beta24、q45のliteral DAGを保存した。
全61とbeta24の原点値が、凍結済みW1認識の出力と厳密に一致することを照合した。
全新norm分母も原点で非零である。
sourceの全Poisson恒等式と座標変換の不変性により、このmapはactual CA族のS61正規化である。
今回134式を再展開するCAS試験は行っていない。

## 2. 旧49座標へは追加のgaugeが必要

旧49表示はq45とcritical beta sliceの4変数(ss,aa,bb,cc)であり、任意のbetaから
単に四つの成分を選ぶmapではない。
Sp4 orbitの7方向を除き、primitive betaのrows[0,4,5,6,7,8,16]を基点の定数へ戻す必要がある。
旧受理済みのorbit minorは-80であり、この7方向はunit formal gaugeを与える。
gauge後に限り

```
ss=beta11/2, aa=beta2, bb=-beta14/2, cc=-beta3/2
```

であり、共通selected10 qとremaining35 qがv,zとなる。
さらに、old49はfull Dから17個のsmooth方向を除いたcritical表示である。
現在のtailはphysical G19切断に加えw=1を課すため、差の三つの残余切断も必要である。
旧218式だけをそのままtailと呼ぶと三方向を失う。

## 3. 新しい20×17単位minor

上のforward gauge全部を解く代わりに、旧49から元physical tensorへ戻すinverse-group方法を調べた。
新W1 frameでvalue補方向の六行列

```
E22,E44,E32,E42,E31,E41
```

とperspective四行列E01,E02,E03,E04、旧Sp4 basisのindices[1,3,4,5,6,8,9]を選ぶ。
全17行列を固定GLで元のphysical frameへ共役した。
同じ保存済みE1 tensorでの全150 coefficient Lie responseを正確な有理数で作り、
共通G19 cut rowsとphysical w=row140の計20行を抜き出した。
この **20×17行列のrankは17、保存した17次minorは32** である。
これは新しいsource jetではなく、同じ固定pointでの明示的なgroup方向の線形係数計算。

従って次の具体的constructive gateが得られた。

1. 旧49の共通polynomial embeddingをS61の61座標へ戻し、150係数を再構成する。
2. 固定GL逆を適用し、上の17方向を接線とする有理group sectionを作用させる。
3. 選択した17個のphysical cut条件を、unit Jacobian32の17-variable Hensel系で解く。
4. 残る三つのcut条件を旧49上のactual関数として保持する。

このliteral group/cut入力を構成すれば、旧218式の商にその三条件を加えたtailモデルと、
physical R0九関数のmapを同時に得られる。old nil Fとその線形annihilatorをそこへ運ぶことも
初めて正当化される。今回これら17未知数の解、残余三関数、nil pullbackはまだ構成していない。

## 4. 次元の取り違えを避ける

旧W1のfull Dはcone dimension28、depth21であり、nilはregular21次元support上の平方零module。
whole ringがregular21次元という意味ではない。
actual family bindingと同時G19 productが完了した場合、G19を除いたW1近傍はdimension9、depth2。
さらにgraded identity
`D_perp[w^-1]=(H_A/Jinf)[w,w^-1]`
から、coefficient tail H_A/Jinf自体はdimension8、depth1となる。
nil tailのsupportは1次元だが、被約N tailは8次元である。
whole tailを`k[[t]][epsilon]/epsilon²`と置き換えない。
この次元計算だけからJinfのpresentationやT3の非零性を得たとは主張しない。

## 保存・実測

`code/ambient_B18_W1_tail_coordinate_preflight_20261104_v63.py`。
`reports/AMBIENT_B18_W1_TAIL_COORDINATE_PREFLIGHT_20261104_v63/tail-coordinate-preflight-DAG.json`に
全53122 scalar node、全150/61/beta24/q45 pairs、valueframeとinverse、全norm分母、
group17の元frame行列と全150×17 response、20×17 matrixとminorを保存した。
一回0.724080秒、最大RSS122620KiB、証明書928427bytes。
GB、source再生成、source jet、追加点、p逆元は無い。
