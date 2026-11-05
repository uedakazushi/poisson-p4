# 残るlower carry：四pair・120 gauge係数のsource preflight

2026-11-05。新しいsource/gauge係数はこのpacketでは評価しない。
実six-adjointがraw X2へ要求する15二次monomialのうち、accepted proper base20に無いものは
`z6*z0,z6*z4,z6*z5,z6*rA`だけである。
base20に既知のfirst z6とこれら四pairを加えると25次元の閉じたmonomial quotientになる。
各monomialの物理source二次値は既存二族own57/source2とfull9 F9 adapterの合成で得られ、追加のphysical implicit solveは要らない。

canonical gaugeのfirst z6も保存済み30 first9から埋める。未知は四pairにおける30係数、計120である。
保存済みJ30 inverseに四RHSを掛ける有限線形補正で得る。補正後の全30条件、base20への射影、実proper κ/G8の同じauxiliary-etaA frameを検証するのが次jobの達成条件となる。
これはgauge係数の新規値なので、parentの別GO前には実行しない。

third sectorのcanonical m1 last-rowにはVgammaのpure-r三次が見えるが、実値は不要である。
`d1*m*Xp=0`をz=0へ制限すると`r^n*(mXp)_last=0`なので、domain QQ[[r]]上でlast成分は零。
first-r係数は`m1_last,r*Xp0=-m0_last*Xp1_r`。
全661で既証の`X1u=Xp0*a1`を用いれば、必要なlast-row contractionは18個の保存済みfirst係数だけに戻る。
従ってこの補完はsource3やproper3を要求しない。

出力は明示したsource/function namespacesと集約slot上界を保存する。正確なHF2や残る四行の値はまだ主張しない。
