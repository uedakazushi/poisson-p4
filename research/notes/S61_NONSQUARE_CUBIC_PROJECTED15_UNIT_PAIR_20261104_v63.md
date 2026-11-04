# genuine二行のheavy pair unit：global cubic pivot とDelta34上のfull57

2026-11-04。FULLU=D(adFBCJ)。個別HZ/HW係数のunit gateが不成立だった有限preflightを保持し、
その実係数二つから明示的な小さなBézout identityを検証する。新点、rank計算、GB、
parameter determinant展開、四次式展開は使わない。これはprimary証明で、独立監査待ち。

## 1. 実source二行とliteral unit

元remaining57 sourceの行42,47はnormalized43,49、actual135 cubic20,61である。
その差をoldH座標 (x0,x1,x2,x3,x4,Z,W) に移すと、最初の五H係数は零で、
H*(kZ Z+kW W)だけを持つ。実係数は

    g1=(3e−2)d+20,  g2=(6e+1)d+40,
    h=Cg1−B,
    kZ=−aFB*g1*g2/10^7,   kW=−FB*h/10^7.

全ての因子式を、保存された実source係数とpolynomial multiplicationで照合する。
L0=6e+1、K=600a+126e²−63e−14−L0*C とすると

    2L0*h+K*g2=240C,
    K*C*g1*g2+(2L0*B−K*g2)h=240BC.

従ってdelta=−64000C²、cZ=K*C、cW=a(2L0*B−K*g2)として

    cZ*(delta*kZ)+cW*(delta*kW)
      =eta=192*a*F*B²*C³/125.

etaはFULLUの既知単元である。個々のkZ,kWを単元と仮定せず、未知因子g1,g2,hを反転しない。
二係数の同時零もFULLUでは起こらない。literal QQ polynomial identityによる結論であり、
有限標数やsampleからのpromotionではない。

## 2. 実global cubic56→90の同cokernel recipe

z=delta*kZ、w=delta*kW とする。H係数column (HZ,HW) に

    T2=[[cZ/eta, cW/eta],[-w,z]]

を作用させると、Cdelta=delta*(G61−G20)の二H係数は(1,0)になる。
det(T2)=1、inverse=[[z,−cW/eta],[w,cZ/eta]]を両方向の積で検証する。
元physical seven H columnsからoldX seven H columnsへの変換は定数可逆行列で、
pure84 columnsを変える必要はない。row coefficient matrixに対してはright T2^tである。

source row47をCdeltaへ置き換え、row42を保持するsource changeのdeterminantはdelta。
他の56行はそのまま。変換後Cdeltaは新Hnu coefficient1、他のH係数0と、元のpure84係数を持つ。
このunit pivotをSchur消去して、**元actual57×91 cubic sourceと同じcokernelの自由56→90**を得る。
これは既存実sourceからの有限recipeで、56×90 coefficient matrixを展開した主張ではない。
quad productsはHを持たないため、この新Hnu方向への係数が零で、元41 product行は変化しない。

## 3. Delta34 chartでのfull cubic57

受理済みactual product41により、degree3 quad商Bはpure43⊕H*light7のrank50 bundle。
ZW制限のfour-quad product19により、その商Vはpure16⊕H*fiveのrank21 bundleである。
したがって0→T29→B50→V21→0はFULLUでexact。
これはdegree3だけの議論で、U78、CI、unmixednessは使わない。

compact relationはprojected genuine16がprimitive line G61−G20をkernelに持つことを示す。
既存の保存restricted36×40 matrixのnonzero34minor Delta34が非零なら、quad19方向を除いた
projected genuine16のrankは15。そのkernelはちょうどこのlineになる。
上のeta-unit functionalはそのlineのfull B50像を非零にし、各fiberでsplitする。
従って他の十五方向と合わせて、fullgenuine16はquad product41から独立であり、
**actual full cubic matrixのrankはFULLU∩D(Delta34)の全点で57**となる。

Delta34は元実polynomial matrixの固定minorで、保存QQ点で非零と既に受理済みのもの。
今回、新点で評価せず、determinant polynomialも展開しない。
Delta34をFULLUの全域単元へ加えない。Delta34=0のlocusでもfull57が成り立つ可能性はあるが、
本packetではそのlocusを処理しない。

global cubic unit pivotはFULLU全体で有効だが、fullrank57の結論は追加Delta34 chart上である。
全FULLU full57、degree4 full239、generic-beta exhaustion、kernel方向の追加は主張しない。
既知uniform kernel38、残り42は不変。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_cubic_projected15_unit_pair_20261104_v63.py --output work/marked-square-v63/cubic-projected15-unitpair-REPLAY
```
