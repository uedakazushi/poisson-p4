# D0一次関係の未知商微分と、source4の限定preflight

2026-11-04。新しいsource四次計算を実行する前に、actual D0の全640関係が
どの未知first derivativeへ依存するかだけを調べる。
T0のsocle、depth、未知微分の値はこのpacketでは決めない。

原点で零である四つの旧CB g商と四つの新CA H=Bbar/gAを、それぞれ独立な
正order functionとして保持した。既存source2から正確に決まる全coordinate商のfirstjetは
使うが、これら八functionのfirstjetは零へ置かない。
旧CB関数は(r_B,s_B)、新CA関数は(r_A,s_A)の二変数functionなので、
最終的な未知derivativeは各二つ、計16個である。

全640×9 entryを一次まで評価し、既知部分と八つの未知functionの係数を保存する。
未知感度の無い列だけの固定spanと、全既知部分・全許容未知derivative効果を加えた
拡大spanを区別する。拡大spanのrankをactual rankと呼ばない。
固定spanに未知効果が入る場合だけ、その部分は未知値に依存しないと判定できる。

source4が必要になった場合の入力も読み取りだけで確認する。
pure CA57方程式、F9、f63と四つのzero physical outputのancestorを数え、
未知数が実際に保存済みCA57内だけであり、全三次source係数が既存joint105内にあることを確認する。
quartic correctionはdegree4の五monomialに対して57×5係数だけ。
全monomial数は二変数degree≤4の15個である。

H(0)=0とBbar=gA*Hから

```
Bbar_4=gA_3*H_1
```

が正確に従う。gA_3のr³係数3/8とs³係数2/3を使う二つのendpoint係数のminorは1/4。
従ってH1の二係数を求め、残る三つのquartic coefficientでも全恒等式を検証できる。
gA4はH0=0のため不要である。
このpreflightは五つのquartic coefficient自体をまだ計算していない。

旧CB商のfirst derivativeがC0表示で相殺するかも同時に保存する。
H9表示で旧CB感度が消えたことを、gA*tauを独立generatorとするC0表示へ無断で移さない。
両族の未知微分が残る場合は、CA source4だけでは全一次行列を決定したことにならない。

Code `code/ambient_B18_T0_firstjet_sensitivity_preflight_20261104_v63.py`。
Report/folder `AMBIENT_B18_T0_FIRSTJET_SENSITIVITY_PREFLIGHT_20261104_v63`。
source solve、新しいjet、GB、point、primeは無い。
