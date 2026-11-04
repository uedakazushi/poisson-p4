# actual T_nil の完全 minimal-target 5×849作用素

2026-11-04。係数環はS=R0/(ell)=k[[zA0,...,zA6,rA]]。
新しいsource/gauge solve、jet、旧r/g函数の再除算は行わない。
独立受理済みrank24/18/81の定数minorを全次数の単元minorへ持ち上げ、完全128×972表示を
同じcokernelの5×849表示へ正確に縮約する。係数級数をscalar展開した表示ではなく、
全actual sourceを固定した有限typed matrix operatorである。

まずEをm0の選択24列とし、n=E(mE)^-1を用いる。mn=I、Vgamma n=d1から既存の
完全比較証明がそのまま成立する。nはliteral E座標にのみsupportedなので、24section単位対を
消すと、上20行が元Xの補集合座標、下84行が元の形式函数lift hとなる104×948表示を得る。
元X0の18選択列はその上20座標でも独立であり、互換な18×18 minorを固定する。
これを消すと上2行と84wedge行、sourceは残る762 X列と全168 d3列となる。
次に固定W0の81行・81列minorを消し、5×849を得る。

各段階はU単元blockに対するD−C U^-1 Bである。target projectionと全source列の経路も
保存する。元780 X関係と全168 d3関係は全て残留列または単元消去列として記録される。
単に定数rankだけを使って関係を落としたのではない。全残留entryの原点値は零で、
五targetはminimalである。これもnilidealの全生成元を与えるという意味ではない。

一回のtyped builderは1.65秒・約25MB、17,716 typed/QQ slots・135KBで完了した。
残るtargetは元free44 source42,43、およびsector0のwedge(0,7),(4,7),(5,7)である。

## d3だけから決まる一次関係

d3の原点値は零で、最初の42単位対の上blockにd3は現れない。最後の81単位対の消去後、
d3の一次値は上2targetで零、残る3wedge targetでは固定constant quotient projectionと
literal Koszul d3一次部分の積である。全168列のこの値とrank証明を保存する。
これは真の一次関係部分空間であり、完全な一次関係空間とは主張しない。
この既知部分のrankは15で、単独ではHF1≤40−15=25を与える。
別の総長上界から得る制限を、未計算の一次関係として付け加えない。

## 全一次層に必要な未計算データ

X0 c=0というconstant residualでは
(mXc)_2=m0 X2 c+m1 X1 cとなる。m2 X0 cは消えるが、m1 X1 cは自動的には消えない。
さらに18pivotの一次source消去係数にはh_pivot,0が掛かり、それもm1 X0pivotを含み得る。
従って旧計算のm1不要という定数rank専用の性質を、一次関係へ延長しない。

m1はordered-prefixのactual Vgamma二次係数を要求し、第三sectorのr二乗除法のr方向には
pure r三乗係数も必要となる。既存C44/B11はr二乗とz0r,z1r,z4r,z5rを保持するが、
全z_i z_j・全z_i rは保存していない。X2はactual adapter、旧B/P、およびradical native35の
同じ二次係数を必要とする。ResAT2*Bold0がconstant residualで消えるかも、元Bold0を
用いて明示的に検査し、結果をprofileに保存する。未知係数を零とはしない。
実際、九つのconstant residualが非零Bold0を保持した。そのspanはrank4、supportは
旧48行の32,36,40,44のみで、physical zero-graph rows64,68,73,74に対応する。
従ってこの項で必要なadapter二次値はこの四列へ限定できるが、全て消えるとは主張できない。

このpacketは上の全次数作用素と既知d3一次部分に限る。全一次Loewy層、総長、
元nilidealの生成元・三乗・正確な指数は決定しない。
