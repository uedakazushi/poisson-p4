# μ5/r1全marked quartic：五有限候補と既知二枝、A3境界control

2026-11-07、PRIMARY E。producer手証明、独立確認前。
数値入口によりmu(F)=5ではr1、全D=HFがquasihomogeneousでなければならない。
全五profile A5,D5,A4+A1,A3+A2,2A2+A1のmarked geometryをここで尽くす。
有限候補の全Saito/cofactor/磁気sourceは次の一括検証に残す。

## 1. A5及びA4+A1

全oneparameter familyからの根比較は `PLANE_HYPERPLANE_QUARTIC_MU5_A5_A4A1_MARKED.md` に証明した。
残る全marked orbitは次の三つ。

```
A5 smooth hyperflex:
 F=(hw-y²-2h²)²-16h²y²-64h4, H=w+6h.

A4+A1 node inflexion:
 F=(hw-y²)²-4hy(h-y)², H=h-w.

A4+A1 smooth hyperflex:
 F=(hw-y²)²-16hy(h+y)², H=h+6y+w.
```

A5 common tangent及びA4 tangentの除外はweighted-branch intersectionの不整合による。
他のline又はglobal modulusの補集合は元noteの全binary quartic比較で保持されている。

## 2. D5の全family：既済cusp tangentと一smooth hyperflex

intrinsic cusp tangentをh=0、transverse smooth branch tangentをy=0に取る。
全integral D5 quarticは、w-shearと非零diagonalsで

```
F=h²yw+h4+d h y³+y4, d全C
```

へ移る。h4,y4の係数が零なら各々y,hがFを割るので両者非零であり、正規化で落とさない。
normalizationは

```
[h:y:w]=[s³t:s²t²:-s4-d st³-t4].
```

cusp/smooth二branchの逆像はs=0,t=0、両者が同じD5点へ写る。
Hがその点を通ってintersection四を持つならh又はyの二tangentだけ。
H=hは `PLANE_HYPERPLANE_QUARTIC_D5_TANGENT_SOURCE.md` で全familyのtau12条件d=0と証明済み。
H=yでは二つのsmooth branchがともにA2 cuspへtransverseで、互いのintersection二となる。
A2 weights(2,3)のsmooth weighted homogeneous branchは二axisだけ。
cuspへのintersection二のaxisは一つしか無いので、二つの別smooth branchは同じweightsでhomogeneousになれない。
従ってDはquasihomogeneousでなく、この全caseは除外される。

HがD5点を避けるならsmooth hyperflex。
general line pullbackは `-C r4+A r3+B r2-C d r-C`。
これが `-C(r-u)^4` に等しいことからu4=1,d=−4u³、A=4Cu,B=−6Cu²。
diagonalの四乗根作用で四caseは同一marked orbit、rational representativeは

```
F=h²yw+h4-4hy³+y4, H=4h-6y+w.
```

F|H=(h-y)^4で、hyperflex点はsmooth。
他の唯一D5とhyperflexのA7からtotal tau=5+7=12。
この一族は既済D5+tangent canonical点と違うmarked branchであり、そのD101/Qzeroを再計算するものではない。

## 3. A3+A2：common tangentだけ、smooth hyperflexは無い

tacnode二逆像をs=0,t=0、cusp逆像をs=tへ取る。
common tangent hのpullbackはs²t²、tacnodeとcuspを結ぶline yのpullbackはst(s-t)²。
第三sectionからh,y成分を除き、cuspでのfirst derivative零条件を課すと全familyは

```
[h:y:w]=[s²t²:st(s-t)²:k s4+2(1-k)s³t+t4],
k!=0,1.
```

k0はbase-point-free/integralityに反し、k1ではtacnodeの二curvatureが一致してA3でなくなる。
全kをこの段階で保持する。
general lineのpullbackは

```
C*k*r4+(B+2C(1-k))*r3+(A-2B)*r2+B*r+C.
```

smooth hyperflexならC*k(r-u)^4で、constantからk*u4=1。
r/r3の比較から1-k=2k(u³-u)、従って

```
(u-1)³(u+1)=0.
```

u1はcusp自身、u−1はk1となり対象外。smooth hyperflexは全familyで存在しない。
A2 tangent contact三は四に足りず、tacnode点で四を持つlineはcommon tangentだけである。

general kは既知new103のF3(a) familyと、
`k=(1-2sqrt(a))/(1+2sqrt(a))` 及びnormalization上の可逆Mobius/base scaleで一致する。
k=−1はa∞に当たるため、既済finite-a source式で割ってはいけない。
その正確なquarticは

```
F=h²(w-2y-4h)²-y³(y+4h), H=h.
```

次のSaito batchではこの一caseを境界controlとして保持する。
general familyのgradient rank十六と、この点のrank十六及びSaito determinant非零が確認できれば、
unitな16minorで全quadratic kernelを正則なrank二bundleとしてliftできる。
任意のspecial anchorframeと全六magがその実族の端点となるのでnew103への閉性を適用できる。
このrank/全source liftの確認前に、divisor閉包だけからspecial source帰属を逆推論しない。

## 4. 2A2+A1：一smooth hyperflex、node inflexionは無い

二cusp inverseをs=0,t=0に取る。joining line hと二cusp tangent y,wをbaseに選ぶと、
全familyは非零scalesを除き

```
[h:y:w]=[s²t²:s³(s+t):t³(s+lambda*t)],
lambda!=0,1/4.
```

lambda0はbase-point/integralityの補集合。
node二逆像u,vはu+v=−1,uv=lambdaであり、lambda1/4では第三cuspへ変わるので別の3A2枝となる。
全node familyを先に代表lambdaへ固定していない。
general lineのpullbackは

```
B*r4+B*r3+A*r2+C*r+C*lambda.
```

nodeだけでintersection四ならroot partitionは3u+v又はu+3v。
r3 coefficient比較とu+v=−1からu0又はv0を強制し、lambda!=0に反する。
二cusp tangentのintersectionは三なのでr1にはならない。

smooth hyperflexならB(r-u)^4。r3からu=−1/4、linear/constantからC=B/16,lambda=1/16、
r2からA=3B/8となる。
従って唯一のwholemarked orbitのrational representativeは

```
F=256y²w²-288h²yw-256h³y-16h³w-15h4,
H=6h+16y+w.
```

F|H=(3h+16y)^4で、hyperflex点はnode/cuspを避けるsmooth点。
元の2A2+A1とA7からtotal tau=2+2+1+7=12である。

## 5. 次の有限batchの正確な全範囲

新たなsource検証を要するfinite候補は、§1の三つ、D5 smoothhyperflex、2A2+A1 smoothhyperflexの五つ。
A3+A2 k−1は別control。既済D5+tangent ambient点は再計算しない。
全18quadratic syzygy・Saito/cofactor・全六magの一batchを準備する。
どのcandidateもtau12だけでPoisson sourceや成分と認めず、全source条件を確認する。
この幾何listと既済mu6listによって全九profileのmarked geometryは尽くされるが、
全E帰属又はcomponent数の尽くしは、source batch・actual認識・独立reviewの後まで未完了である。
