# reduced比較余核のgeneric DVR rank2：独立幾何監査

`Cprime=coker(D0→E4_red)` は同じactual `C=coker(D0→E4)` の商である。
従って既受理のDVR-finitenessと4曲線のsupport包含が移り、新しい支持枝は加わらない。
各曲線のell order1・DVR residue degree1も、既存の同じsource embeddingの独立監査を使う。

actual4写像をそれぞれ被約targetへ降ろす。CA-criticalではDredとEredが同じ
`G[[l7,v]]`、actual写像は恒等なのでcoker0。Paでは`O_N→O_N⊕O_SP`が制限graphなので
coker O_SP、generic長さ1。CB-criticalはactual CAの外なので0。
非対称Qbではaccepted実写像
`B/(fJZ)→B/(f)⊕B/JZ`、`JZ=(z7,f²,rho²)` を使う。
被約targetは `B/(f)⊕B/(z7,f,rho)`、第二写像は第一の商なのでcoker長さ1。
Qbのnil長4や元coker長2と、被約coker長1を混同しない。

graphのcokernelは `(a,b)→b−pi(a)` で直接同定する。これらは抽象模型を後から
alignした写像でなく、各whole proper normalization source diagramに束縛された実写像。

zero比較への輸送にも条件が要る。generic零次数conductorでは相手conductorとpositive Jが
単元になり、full Cをwが殺す。さらにfull Cの生成元tauA,tauB,tauA*tauBはweight0なので、
full Cはdegree0に集中する。w-annihilationだけではこの集中性は推論できない。
full reduced-target cokerはそのgraded quotientで同様に集中し、非負gradingから
`(Efull_red)_0=(E4)_red`、degree0のexactnessから同じCprimeが得られる。
actual G19群積はnormalizationと比較写像を同時に保つ平坦な正則因子であり、
任意の非平坦なw断面を選んだ議論ではない。

従って寄与 `0+1+0+1` とfield degree1よりrank_D Cprime=2。
独立値監査で `length(Cprime/ell)=4` だからDVR最小生成数4、torsion最小生成数2。
非標準的に `D²⊕D/(ell^m1)⊕D/(ell^m2)`、mi≥1であり、指数・torsion長は未知。
有限商のr-chainはDVR parameter ellとは異なり、指数を決めない。
元nilideal生成数、K/JK、cube、全域支持の尽くしは主張しない。

再現：`.venv/bin/python code/ambient_B18_reduced_comparison_generic_rank2_review_20261104_v63.py`
