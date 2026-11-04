# conductor transfer 512 actual cycle：独立 source 監査

対象は元B18零次数交叉D0から正規化比較E4への実比較。producerをimport・再実行せず、新しいsource jet、scalar512展開、kernel solveをしない。25,088 typed nodeの全型・外部namespace・hashを照合し、全512本のC0 free9 representative、九つのC係数、完全44cycleを別のsemantic fingerprintで再構成する。

## 全conductorとsource frame

CAはactual C0のβ=(zA7,gA)を使う。CBはphysical rows[18,19,28,29,60,67,72,63]のraw八係数を保つ。F_jのτA係数が正、τB係数が負であることを保存されたscalar-function polynomial DAGで確認する。最後のraw63係数はgBと最初七係数のunit-unipotent変換であり、det=1、除外する未知因子はない。両側各四つのzero-graph membershipも既監査sourceに束縛する。

C/R0の基底は(1,τA,τB,τAτB)。τA²のsourceはactualC0 node20086で、pA²=sA*pA+tAからτA=pA−sA/2としてνA=tA+sA²/4となるliteral DAGを確認する。C0のodd八列がβ_iτAであることも全列で再確認する。nil二生成元n1=dh,n2=dclはactual native radical groups5,7、全四基底倍はQrad列4g+γに一致する。

actual12→9adapterの向きは列生成元F=A(dR9)。転置を挿入しない。p=nαeγ=(p0,pA,pB,pAB)、F_j=δ+bAτA−bτBのとき、独立な十変数疎多項式環で

`U(E,H)=β*b*p+β*F_j*(pB+τA*pAB)`

の四成分を確認する。E=β(bp0+δpB+νA*bA*pAB)、H=bpA+bA*pB+δpAB。Qradの実source係数は九つのC係数 `β*A_jg*(pB+τA*pAB)+[g=gα]β*b*eγ`。これを四成分へrestrictした36係数とHのodd八係数が、保存された各44cycleに全て一致する。even成分の回収も同じ四成分恒等式による。既受理のcomplete radical comparisonにより、これは実nil classを表す。

## transfer algebraと生成範囲

各conductor pair(a,b)のφ:E4→D0はD0-linearで、ψφ=f id、φψ=f id、f=ab。同じfの別分解でφが一致するとは仮定しない。一般に `φ_f(x)φ_g(y)=gφ_f(xy)=fφ_g(xy)`、特に同じφで `φ(x)φ(y)=fφ(xy)`。

NE³=0なので各φ(x)は三乗零、従ってW=Σφ(NE)はN=Nil(D0)の部分ideal。さらにψ(N)⊂NEとD0-linearityから `W*N²=0`。`ab*n=φ_(a,b)(ψ(n))` より `(IAIB)N⊂W`。従って三重積写像は正確にSym³_D0(N/W)を経由し、IAIBはN/WおよびN³をannihilateする。cubeそのものが零とは言わない。

NEは二ideal生成元と四R0-module基底の積で張られ、conductor八×八とtransferの線形性を合わせ、512本は上記W全体を生成する。N全体の生成元、最小生成系、非零性を求めたとはしない。

## K/JKの別gate

Jplus=IA+IBとJtimes=IAIBを区別する。既受理JK²=0からpairing K×K→K²はK/JplusKへ因子化。D=k[[ell]]上のrank4とμD K≤13、ell²∈Jplusよりdim K/JplusK≤26。13生成元の対積91個とell²annihilationからlength K²≤182。R0/Jplus²の長さ31も7個のu一次×(1,r,s)の21個とpure(r,s)次数≤3の10個から得る。これらは全H0やN³の長さ上界ではない。

free comparisonをJplusでbase changeしたhomologyとK/JplusKの間にはTor₂(C,R0/Jplus)のconnecting imageが残る。完全なsource cycle liftまたはhonest K-presentationなしに両者を等置できない。Artin–Rees strictness、nilcube、正確なnil指数は未決のまま保持する。
