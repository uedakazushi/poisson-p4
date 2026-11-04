# B18 T3の次の有限入力：完全882 syzygyと小さいD4解像度への経路

2026-11-04。既受理のhonest9単位消去を再実行せず、有限module計算の方法だけを
前進させる。新しいsource jet、点、素数、GB、module membership計算は無い。
**T3の零・非零と、元B18全Dの深さ19/20は未決定である。**

## 1. 今回の出発点と、特殊化で失われるもの

実際の共通係数環をR=R0=QQ[[z1,...,z7,r,s]]と書く。これはactual CAの
F9 unit座標であり、各成分のtraceを同一視しない。
`AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63_REVIEW`は、10/3の単位で
正確に一行消去し、

```
M=D3=coker(A:R^2069 -> R^9),
N=D4=coker((A,P):R^2081 -> R^9),
T3=image(P:R^12 -> M)=ker(M -> N)
```

を独立受理した。これは二行trace ambientへの像だけではなく、元の全36 Koszul
関係を含む正確なmodule表示である。mu(M)=mu(N)=9。

両側actual zero odd idealの和について既受理の正確な等式は

```
R/(I0_CA+I0_CB)=QQ[[r,s]]/(r,s)^2=:A3,
(I0_CA+I0_CB)*T3=0.
```

従ってT3は12生成、長さ36以下である。しかしT3をMの中へ埋め込んだまま
z7=0やm_R^2=0にする操作は単射を保たない。
例 `z/(z^2) subset R/(z^2)` の像はmod zで零になる。
HONEST9_LOW2は正確な表示からの負試験として正しいが、top12が二次まで零だった
結果は、この構造的な問題を解かない。今回はそのproducerを独立再実行しない。

## 2. CB99の全first syzygyを明示的に作れる

別のactual CB regular coordinatesを用いて

```
b=(z1,...,z7,g,wr,ws),  g=r*A+s*B != 0
```

とする。ここでA,Bは既存の正則cofactor、gはw,z7に依存しない。
coreのHilbert--Burch matrixは

```
[[ w, 0 ],
 [-A, s ],
 [-B,-r ]].
```

core ideal `(g,wr,ws)` のgradeは2であり、七つのzはその上でregular。
従ってKoszul(z7)とcore HB解像度のtensor productは正確で、ideal(b)の最初の
自由加群の大きさは10,44,112となる。44本のlinear relationをS_alphaと書く。
CB algebraは、係数nu=nu_CBと十odd変数Z_iについて

```
S_alpha=sum_i s_alpha,i*Z_i,             alpha=1,...,44,
P_ij=Z_i*Z_j-nu*b_i*b_j,                0<=i<=j<10
```

の全99関係で定義される。この99関係の**全first syzygy**には、次の882本という
具体的な生成系がある。

1. Koszul(z7) tensor HBの次の微分から112本。
2. 全alpha,kについて `Z_k*S_alpha-sum_i s_alpha,i*P_ik=0`、440本。
3. monomial ideal `(Z)^2` の三重overlap、330本。

第三群の補正は、`K_ik=b_i Z_k-b_k Z_i`を44本へ戻して作る。
regular/crossのKは既存の一列そのもの。coreの三つは

```
K_g,wr = -r*S_H0+B*S_H1,
K_g,ws = -s*S_H0-A*S_H1,
K_wr,ws = -w*S_H1.
```

例えば

```
Z_k*P_ij-Z_i*P_jk+nu*b_j*K_ik=0.
```

独立の疎整数多項式演算で、全882本と45個のKの戻しを全係数検証した。
g=rA+sB以外の特殊恒等式を使わないので、actual coefficient functionsへ代入できる。
相対w次数は、99関係が64/32/3、882 syzygyが448/368/64/2である。

### 生成性の証明

Z次数を使う。任意の99関係間のsyzygyについて、第二群でS_alphaの係数に含まれる
正Z次数をP項へ移す。残るP項の最高次数が相殺する場合、その先頭部分はmonomial ideal
`(Z)^2`のsyzygyであり、三重overlapで生成される。
第三群を引くと、そのZ最高次数は下がり、補正されたS項も第二群で再び処理できる。
最後はR係数だけのS間のsyzygyとなり、第一群112本で尽くされる。
これは有限次数降下で、flatな特殊化や数値rankを仮定しない。
従って全z7-syzygyも失われない。

87 physical graphはmonicな追加target座標の消去である。
CB99の解像度とその87変数のKoszul complexをtensorすれば、全actual CB idealの
解像度になる。first free ranksは1,186,13236。
13236列を展開する必要はなく、882 coreとKoszul87をfactoredなまま保存できる。
ただし、これが現行9行presentationより小さいという主張はしない。

## 3. 正確なw-Bocksteinと、残る計算

共通full target環上の上記CB解像度をF_CBとし、actual CA algebraへtensorして

```
C=F_CB tensor A_CA
```

とする。H0(C)はactual Dである。A_CAでwはnonzerodivisorなので、Cの各項でもwは
nonzerodivisorであり、complexの完全列

```
0 -> C(-1) --w--> C -> C/w -> 0
```

がある。そのdegree4の長完全列から正確に

```
T3 = coker[ H1(C)_4 -> H1(C/w)_4 ].
```

接続写像は、`d1(c)=w*q`であるcycle c mod wを `[q] in D3`へ送る。
必要なcomplex differentialはd1,d2までであり、上でその全component入力を与えた。
coefficient moduleも、正次数は `R+m9*tau_CA`、零次数は
`R+(z7,g_CA)*tau_CA`なので、それぞれKoszul9、Koszul8の既知の完全解像度で
自由加群へ持ち上げられる。

**H1(C)を零と置いてはいけない。** これは一般にderived intersectionのTorである。
CB99のsyzygyをCAへ代入しただけで、actual Jinfの全syzygyが得られるわけではない。
このBocksteinは、計算すべき二つのhomologyと比較写像を明示する。
今回、それらのkernelやimageを計算したとは主張しない。
大きいfactored complexを無条件に展開して次のjobにすることも勧めない。

## 4. 新しいliteral tailから、より小さい解像度を作る手順

新しい `AMBIENT_B18_W1_TAIL_BINDING_20261104_v63` は、同じphysical R上で

```
0 -> L -> N=D4 -> N_tail -> 0,
L=(R/P_curve)*F_old, R/P_curve=QQ[[T]]
```

を指定する。これは今回独立再監査した入力ではなく、rootが受理したliteral bindingを
使う候補経路である。P_curveはregular8生成、uniformizerはphysical row117、
21-unit nilcurve系のdeterminantは-80である。全D4を一変数二重環へ置換しない。

次の**actual同時座標gate**が満たされれば、小さいN_tail解像度を明示できる。

```
R = R8[[v]],   R8=QQ[[q1,...,q8]],
N_tail = R8 + m8*tau,   tau^2=nu in R8,
v=q1*tau.
```

必要なのは、現在の九physical base functionsのうち八even functionがR8を与えること、
残りがeven部を引いた後b*tauとなり、bがR8のregular coordinateへ取れること。
individual Ac7型の一致や、一次rankだけからこの全同時座標を仮定しない。
actual source/gauge/cut mapによるunit IFTの確認が必要であり、これは別packet待ちである。

### 4.1 N_tailの明示的cone

R8上ではN_tailはgenerator `E0=1, Ei=qi*tau (1<=i<=8)`と全28 Koszul関係を持つ。
m8のKoszul解像度へR8のfree even summandを加え、その上でvの作用をchain map Tで

```
T0(E0)=E1,
T0(Ei)=nu*q1*qi*E0,
Tj=0  (j>=1)
```

と持ち上げる。T0は28 Koszul関係を文字どおり零へ送るので、これは全chain mapである。
Rへ係数拡大し、`v-T`のmapping coneを取る。
作用は位相的nilpotentで、module上ではT²=nu*q1²。従ってmonic module divisionにより
coneはN_tailを正確に解像する。

初めのrankは9,37,84。row E1 / shifted E0列のunit -1を消去して

```
F0_N=R^8, F1_N=R^36, F2_N=R^84
```

を得る。独立polynomial checkerはこのcone、unit消去、全d1*d2=0を検証した。
これは条件付きの普遍的計算であり、actual N_tail gateの計算ではない。

### 4.2 Nilcurve拡張のhorseshoe

Lはregular8 ideal P_curveのKoszul解像度を持ち、初めのrankは1,8,28。
上のN解像度とactual extensionからhorseshoeを作ると、D4のfree ranksは

```
F0_D=R^9, F1_D=R^44, F2_D=R^112.
```

これは解像度の上限であり、higher differentialの最小性は主張しない。
構成すべき新しい係数も具体的である。

* Nの36関係をactual D4の九generatorへ持ち上げる。残差は一意に
  `e_i(T)*F_old`である。
* Rからnilcurveへの写像は全射なので、e_i(T)をT=row117の関数としてRへ持ち上げる。
* d2_Nとの積がLで零になる恒等式を、P_curveのKoszul分解で持ち上げる。

これにより、符号を一貫させて

```
d1_D = [[d1_L, e], [0,d1_N]],
d2_D = [[d2_L, h], [0,d2_N]],
d1_L*h+e*d2_N=0
```

という明示的matrixを得る。未知のSyzJinfを呼び出す処理ではなく、regular sequenceの
既知Koszul divisionと一変数extension係数の計算である。
actual F_oldのH_A-pair liftとこれらe,hは、まだ構成されていない。

## 5. この小さい解像度ができた後のT3の有限presentation

ここは方法を単にSyzJinf criterionへ戻さないための具体的な出口である。
mu(D4)=9であるため、actual九generator間のcomparison mapを原点で可逆に選べる。
小さいd1_Dの44列は、共通R^9からD4への全kernelを生成する。
各old D3 relation A_jをこの44列へ持ち上げたvectorをa_j in R^44とする。
すると**正確に**

```
T3 = coker( [d2_D | a_1 | ... | a_2069] : R^(112+2069) -> R^44 ).
```

証明は、d1_Dの像がker(R^9->D4)、その全syzygyがim d2_Dであり、
T3はこのkernelをold D3 relationsで割ったmoduleであることによる。
comparison liftsの選択差はim d2_Dに入り、同じcokernelを与える。

liftの実行にも具体的順序がある。まずN側へ射影し、上のmonic v-divisionと
Koszul8でN関係へ分解する。そのD4 liftを引いた残りはnilcurve側なので、
P_curveのKoszul divisionで分解する。旧regular r/gを零 locusで0/0にする必要は無い。
ただしexact liftを得る前にそのoriginやjetだけで代用してはいけない。

この**完全presentationを得た後に限り**A3へbase changeできる。
T3自体がA3でannihilateされるのでright exactnessにより同じT3となり、
44行・2181列の全entryはQQ基底1,r,sの三係数へ落ちる。
QQ-vector-spaceで最大132行の有限linear algebraとなる。
既知top12生成を使えば最終module長は36以下である。
これはambient imageをA3へ下げる操作とは異なり、全112 syzygyを先に保持する。

従って次の有用なjobは、高いsource jetではなく、**actual Nの同時R8 gateと
F_oldのactual H_A lift、36個のextension residueを同じnilcurveへ結ぶ入力**である。
現時点ではそれらは未完成で、この小さいpresentationをactual T3へ適用していない。

## 6. 現在のcoefficient modelを小さく持つ際の注意

498k scalar nodesは498kの自由変数ではない。honest10のdistinct variable namesは
1473個で、内訳はphysical base9、二つのoriginal48-unit source、24個の57-unit pointed
inverseである。未使用aliasを含むvariable node数1530とは区別する。
regular coordinate quotientは1576個、非coordinate g quotientは4個。

coordinate quotientは、二つのunit inverse間のdivided-difference matrixを使えば
nonunit inverse無しに表せる。F(U,x)-F(U',x')=A(U-U')+c*coordinate、det A(0)がunitなので
`(U-U')/coordinate=-A^(-1)c`。output差にも同じtelescopingを適用する。
これはunit系の有限matrix recipeであり、新しいsource変数の任意branch選択ではない。
今回は1576個の全recipeを再構築してはいない。

四つのg quotientについては、既証regularityだけを使う。equation `g*h=b`だけを
足すとg=0上の偽成分を作り得るので、無飽和のpolynomial-Henselモデルと同一視しない。
今回gの逆元やbranchを導入しない。この点でも、小さいactual D4 extension routeが優先される。

## 保存と範囲

Code: `code/ambient_B18_T3_syzygy_bockstein_preflight_20261104_v63.py`。
Report: `reports/AMBIENT_B18_T3_SYZYGY_BOCKSTEIN_PREFLIGHT_20261104_v63/report.json`。
全99 generic polynomial、全882 sparse syzygy、conditional N coneを同じfolderへ保存する。
疎整数演算と入力supportの読取りだけ。全actual CB99の生成性は既受理入力を再利用し、
今回の新しい882生成性は上のZ次数降下証明による。
T3のactual presentation完成、非零性、depth19、depth20のいずれも主張しない。
