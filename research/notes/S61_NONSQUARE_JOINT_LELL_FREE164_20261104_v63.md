# quad117を先に消し、正確なsource補空間を構成する

既存の mixed free183 行列から19本を単に取り除くことはしない。
元319 source で二次の98+5+14 pivotを先に行う。
これで source202=quad57+nonquad145、target122の損失のない presentationになる。
全33genuine quarticはこの段階で保持する。

二次の行操作のみなので、旧15方向の cubic/quartic pivot と新separatorの
source係数は変わらない。旧15を正確に消すとquad57+surplus130。
quad57の22本のkernelを捨てるには**実左逆**を構成する必要がある。

純二次kernelの係数を57座標から復元する。forward行操作が
`g_r<-g_r/p, g_t<-g_t-c_t*g_r`なら、係数の逆輸送は

```
k_r=(k'_r-sum_t c_t*k'_t)/p.
```

14段・5段を逆順に戻し、その後既存98段 restoration を合成して旧174へ戻す。
overlap5でのq0/q5の58本の係数は、既存の正確なNF98 circuitから5列だけを取り
再構成する。genuine cubic/quartic係数の展開はしない。
既存P22をこのrestorationに作用させたものがP57。
K57は変換した旧K22を残る57座標へ制限したもの。
純二次恒等式と117 unit pivotからrestore(K57)=旧K22なので、P57*K57=I22。

Bnew22はquad57から選ばず、surplus130からseparator二座標を避けて選ぶ。
`[A57+Bnew*P57,C108]` は自由165 source。C108に残る既知単元separatorで
constant pair changeと一削除を行うと自由164→122を得る。
この実source basis構成により、旧319と同じcokernelである。
単なるrank差計算やprojective quotientの自由性の仮定ではない。

全復元係数、P57/K57と新kernel165、source/target全index、target NF239→122を
有限scalar/linear DAGで保存する。全164x122を展開しない。
既知uniform方向38・残り42は不変。全点rank122/full239は未解決。
