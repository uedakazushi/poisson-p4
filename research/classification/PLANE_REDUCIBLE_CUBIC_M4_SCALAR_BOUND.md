# reduced cubic plane：m4 squarefree・非特別crossratioの三型G上界

2026-11-07、主義務 G。独立review待ちの今回の証明。
Rがtransverse conic+line、tangent conic+line、triangleの三型で、
Dのplane沿い位数m=4、leading binaryquarticがsquarefree、classical invariant I,Jがともに非零とする。
この全支持のprojective PGL像はそれぞれ次元≤24,23,23。
concurrent型、I=0/J=0、non-squarefree、m2/m3は未完了として保持する。

## 1. 全leading symbolがnormal接続を有限化する

`PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md` のHamiltonian associated-graded証明により
γ4=G4/f は conormal接続にparallel。G4はbase一次係数のbinaryquartic。
相対GL2作用で I,J のweightは2,3なので J²/I³ は不変、従ってconstant。
G4の係数次数から I(G4),J(G4) はbase degree2,3のpolynomial。
UFDにより `I=aL²,J=bL³`。squarefree性よりΔ(G4)=cL6、c≠0。
parallel tensorのdiscriminantはU=S\Rで零にならないため、LはRの一次成分である。

binaryquarticの四rootが相異なりcrossratioがconstantなところでは、三rootを局所holomorphicに
正規化するprojective frameがあり、そのframeの接続が唯一のprojective接続である。
よってG4のdiscriminantが非零なRの他成分上でもtracefree接続は正則に延びる。
元のΩはsimplepoleなので、唯一の残る極Lを持つglobal section H0Ω1(L)⊗sl2となるが、
H0Ω1(1)=0。従ってΩ=0。
ここでnonlog一般枝を捨てていない。この指定leading-symbol枝に限り、極supportの次数から零と証明した。

projective connection零によりG4は `L*constant binaryquartic`。
normal bundle scalar connectionをβとする。
canonical P_f のnormal座標のcube/fはCasimirなので、local frameでβ_can=(1/3)dlogf。
γ4=(L/f)g4 のweightは4−2=2だから、parallel性は
`dlog(L/f)+2β=0`、従って

`β−β_can=−(1/6)dlog(L³/f)`。

これで全normal接続が固定される。m>2なのでnormal-normal quadratic Kは零。
従ってcanonical Euler gauge後のweight0は一つの
`P0=P_f+Euler_normal∧V` に固定される。
Vは上の閉scalar一形式から得るlinear base fieldであり、追加parameterはない。
Rのcomponent-line選択は各型の有限automorphismで入れ替わり、一代表で全てを覆う。

## 2. 全normal weight1/2と三角fiber上界

既済のglobal filtrationにより、全残源はweight1三十係数、weight2九係数だけ。
その完全basisを用い、全350projective Jacobiのweight1と、D3=0の全24係数を構成した。
三型ともそのkernelは八次元。
weight2の式は

`D2(P2)+Jac(P1)=0`

で、未知の全九P2に関するlinear operator D2はP1によらず固定。
三型ともrank6である。従って各任意P1についてP2解は空又はaffine三次元。
higher-weight Jacobiは更なる条件なので、解を増やさない。
全fixed-P0源の次元は `8+3=11` 以下。
未知fiberのgeneric rankを全点へ推論する議論ではなく、同じ固定linear operatorを全P1に用いる。

plane+Rを保つconnected PGL subgroupの次元は
`dimAut(R)+normalGL2(4)+normal→base(6)=d+10`。
dは順に1,2,2。normal接続はRとLの内在的な式から固定され、connected groupはLも保つので、
tensor全体をf係数1へ正規化した全源をこの群が保つ。
したがって `PGL5×X→P H0(Λ²T)` のfiberにこの群が入り、像次元は

`24+11−(d+10)=24,23,23`。

## 3. 有限証拠とconcurrentの未完了

契約は `CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SCALAR_JOB_20261107.json`、
reportは `CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SCALAR_V2_20261107.json`。
独立Fraction全350行、三十・九の全basis、全kernel、全source矩形matrixを保持する。
各P0自身の全Jacobi零、actualD零、scalar接続の全linear係数も検査した。
0.7055秒で全算術PASS。

concurrentではweight1kernel14、weight2rank3/fiber6、粗像上界31。
V2で全G4=L*g4 leading-shape条件と全P1二次polarizationを追加しても、linear P2rankは3のままである。
この結果はreportでfalseを保ち、G≤24を主張していない。
二回の有限linear-bound試行後、この経路は停止した。
次の問題はP1の非零二次compatibilityとD4squarefree条件の構造であり、
minorやjetを増やすだけの探索へ進めない。
旧default/旧m4reportは当時のsourcebytesごと保持する。
