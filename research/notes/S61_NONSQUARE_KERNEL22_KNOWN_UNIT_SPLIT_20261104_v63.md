# Frozen K22のglobal20 splitと残余154×2

2026-11-04。基礎開はFULLU=`D(a*d*F*B*C*J)`。既受理174×22のKだけに、
一回のknown-unit Gaussian前処理を適用した。元quadやquartic matrixは再構成せず、
未知因子を逆にせず、minor/GB/新点の探索はしない。

結果は **20 global unit pivots**。native0.342秒、compile等を含め全体1.746秒。
個別known-unit entryが無くなったところで停止した。P22やglobal282 basisを
完成したとは主張しない。追加のpivot/minor/Bézout jobは実行していない。

## 1. exact20 left split

選択したoriginal source174 rowsをsigma、kernel22 columnsをtauとする。
保存A=`K[sigma,tau]`は20×20。nativeは全係数で **A=L*U** を検査した。
Lはunit lower triangular、Uはupper triangular、その全diagonalは
非零QQとa,C,Jの整数冪だけからなるunitである。一般に許したunit listは
`a,d,F,B,C,J`。unknown polynomialを割らない。

`P20=U^-1*L^-1*pi_sigma : R174→R20`

を、forward/back substitutionによる2044 nodes/73 scalar atomsの有限exact回路で
保存した。pi_sigma以外の154 source columnsの係数は零である。
従って **P20*K[:,tau]=I20**。これは全22列へのleft inverseではない。
source/source-kernelの順序、元A/L/Uの全有理係数、unit pivot logを保持する。
inverseのscalar denominatorを無制限にfieldから選んだものではない。

一般的な方法として、quad coefficient matrix Q6×28にglobal right inverseVが
得られれば、cofactor matrixX6×28を `(XV−(XV)^t)/2` で15 Koszul座標へ送ると
Koszul15のleft splitになる。ただしquadのpointwise rank6だけを固定minorのunitと
取り違えてはいけない。今回はそのVを新たに構成せず、実Kに許可されたunit操作を
適用して20方向を同時にsplitした。

## 2. actual fixed M284

前packetの15 fixed SOURCE deletionsは、元319から同じcokernelを持つM304を与える。
今回の選択20 kernel directionsもsource coordinates sigmaでunit minorを持つ。
従ってそこに対応するquadratic20 rowsをさらに単に削除できる。

`certificate.json` の `retained_actual_source319_indices` が、同じ実像/cokernelを
持つ **固定284×239 submatrix M284** を指定する。新しい全quartic係数は一つも
展開していない。37方向の独立性という前の結果を変更せず、そのうち35方向を
global固定座標で除去できたという進展である。

## 3. 残るactual154×2とcase cover

未選択kernel columnsは

* column12: 元quadricsのKoszul(3,4);
* column21: 新しい22本目

である。K_Iを選択20列、K_Jをこの2列とすると

`K'_J=K_J−K_I*P20*K_J`

はsource sigmaで零になり、残る154 source座標によりexact **154×2** 行列を
与える。保存native Schur stateとread-only frontier readerでその全係数を保持する。
元22本のpointwise independenceと上の可逆column操作から、この2本もFULLU上で
pointwiseに独立である。個別entryのunitが見つからないことは、unitな2×2 minorや
unitなentry combinationが無いことを意味しない。

全ての2-row minorsのopensはFULLUを覆う。それぞれのchartで残り2 source rowsを
消せば元の282 atlasへ戻る。ただし今回その2×2 minorsを計算しておらず、
global P2/P22・単一のglobal282 submatrixは未構成である。

既に計算済みSchur stateの読取りだけによるfrontierは、active106 rows、両列active89、
非零195 entries、7763 numerator/denominator terms、numerator degree最大14、
denominator degree最大3である。全active pair数5565。両列activeな最小term12行から
66個のcandidate pairを**列挙するだけ**で、determinant・factorization・Bézoutを
実行していない。candidate listは不存在一般論ではない。

## 4. 数学的な意味と残るgate

元のrank282 source bundleは今や `coker(R2→R154)⊕R130` とも表せる。
global282 basisは未完成だが、固定M284は既に使える。そのtranspose cokernel D45も
rank45/pd≤1で `Ext¹(D45,R)=actual fullquartic cokernel` を満たす。既存のrank43
双対表示と同じ欠損を扱う。codimension2だけで全射とはならない。

次の有限目標は、この154×2に対する既知unitの2-minorまたはentry/minorの
明示的Bézout splitである。pointwise rank2はcase coverを保証するが、固定unit minor
の存在を保証しない。fullquartic239、actualC122 rank64、残る43 kernel directionsは
なお未解決である。

実行はfrozen K sourceから一回だけ。native checkpointは毎pivot保存し、
35秒探索/全体60秒2GB/300k terms/32MiB上限を指定した。実際の最大Schur termsは
log内の46973、最終9149、capに達していない。検索を自動拡張していない。
