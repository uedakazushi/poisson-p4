# double_A全六磁気源から既知cuspidal24への実一致：独立review

2026-11-07、PRIMARY E_multiple_hyperplane。producerと別作者の全source幾何・独立Fraction算術。**PASS**。
**double_Aの全六磁気源、その全PGL飽和及び閉包は、既知成分
R4_GA9_CUSPIDAL_PLANE_QUARTIC_24に入る。**
新成分は加えず、nil22 imageのE未帰属はstable二源だけとなる。

## 1. 入力と有限契約

固定入力は公開 `NIL22_DOUBLE_A_CUSPIDAL_LIMIT.md` とそのreport、
既済 `CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json` の
`actual_family_power3`、及び元の
`CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json` case0の
**unscaled** `full_six_magnetic_source` である。
parentの全七sourceと既知24成分性、元の全十一源入口は既済の独立証明を明示的に依存とし、今回再証明していない。

全hashと90秒/1536MiBの有限jobは `jobs/nil22_double_a_cuspidal_review.json`。
独立checkerは `code/classification_nil22_double_a_cuspidal_review.py`、結果は
`reports/CLASSIFICATION_NIL22_DOUBLE_A_CUSPIDAL_REVIEW_20261107.json`。
gate PASS後、0.267秒・RSS14,284KiBで全恒等式PASS。
producer、producer helper、SymPy又はCASはimportしない。
既存の別作者Fraction(ε)及び疎多項式算術だけを再利用した。
初回の比較は省略された零bracketと明記した零bracketのdict形の差で停止したため、
全十pairを零補完して比較する形に修正して再gateした。非零係数の数学的discrepancyは無い。

## 2. 全anchorからactual frameを別に解く

parent coordinatesを(h,x,y,z,w)、targetを(h,q1,q2,p1,p2)とする。
parentのplane quarticは

`F_a=(hw+6y²)²−144y³(h+ay)`。

a=1/4では `F=hC3`、`C3=hw²+12wy²−144y³` と独立に展開した。
target cubicを

`C=(384p1²p2+512hp1²+384hp1p2−9hp2²−432h³)/384`

とすると、base変更

`y=−4p1−3h`, `w=−48p1−9p2/2−54h`

のdeterminantは18で、C3=−864Cを全係数で確認した。
この因子形の一致だけで帰属は推論しない。

parent/targetの全四anchorを公開sourceから別々に作り、
`A M_parent B^{-t}=M_target` の全quadratic coefficientによるoverdetermined QQ線形系を独立Gaussian消去で解いた。
全equationを満たす一意な解は

`A=[[-1/448,−5/1344],[11/252,−5/252]]`, `det A=5/24192`。

baseとfiberを合わせたparent→targetのactual行列をこの解から再構成した。

```
G=[[1,0,0,0,0],
   [0,-1/448,0,-5/1344,0],
   [0,11/252,0,-5/252,0],
   [-3/4,0,-1/4,0,0],
   [-4,0,8/3,0,-2/9]],
det G=-5/435456.
```

全inverseとG G^{-1}=Iを再計算した。定数で可逆なPGL5であり、端点へ極を持つframe変更は無い。

## 3. 六磁気の完全な正逆変換

parent全磁気項はR=r0h²+r1hy+r2hw+r3y²+r4yw+r5w²。
target磁気M=m0h²+m1hp1+m2hp2+m3p1²+m4p1p2+m5p2²へのmapを、
全二次monomialのbase変更から独立に計算した。
正確に

`M=(5/24192)R(h,−4p1−3h,−48p1−9p2/2−54h)`。

全6×6行列はproducerの値と一致し、そのdeterminantは

`(5/24192)^6 * 18^4 = 15625/1909592158002053382144 !=0`。

逆行列と全六式r(m)も独立に作り、公開 `actual_parent_r_binding` と一致した。
どの磁気係数も反転せず、全零磁気を含む全C6の同型である。

## 4. 全150族、全350/strict10、正しい実Pfaffian

a=1/4+ε、r=r(m)とし、元のparent sourceの全十bracketへ
実際にG^{-1}座標と全r(m)を代入し、Gでpushforwardした。
この `G_*P_(1/4+ε)(r(m))` は保存された全rational familyの全係数と一致した。
ε=0を全係数へ代入すると、元のunscaled double_A全源と**全150係数で**一致する。
fiber一次括弧もこの計算に含み、anchor又はPfaffianだけの一致で済ませていない。

全mを独立変数としてstrict homogeneous Jacobi十成分及びprojective Euler-Schoutenの全350係数を
別に展開し、全て恒等零と確認した。
正しいEuler-contracted Pfaffianをcurveから直接作ると

`D_ε=5h/(36ε−5) * [hC+(ε/6)(4p1+3h)^4]`

であり、公開式と一致する。
従って端点は `D0=−h²C !=0` で、全磁気値についてrank四の非零projective tensorである。
producer内の途中のdeterminant変換診断を幾何恒等式として使わず、実tensorから直接求めたPfaffianを照合した。

全rational coefficientの分母は192(36ε−5)を割り、ε零で単元であることを個別に確認した。
parent frame分母18a−7も端点で−5/2。
quarticのprimitive unit1−4aは端点で零になるが、実tensorの分母ではない。
generic irreducible-quartic開集合がこの点を含むとは主張しない。

## 5. 全源から既知成分への正しい矢印

任意のtarget mについて、全六r(m)を持つgeneric ε点は既済のwhole七sourceに属し、
従って既知24成分Yに入る。
上のliteral rational familyはε零で正則かつ非零だから、Yの閉性によりその端点もYに入る。
端点自体が元の全六源に正確に等しいため、稠密な単一orbitから磁気零へ追加推論する必要も無い。
YはPGL不変かつclosedなので、全sourceの飽和と閉包もYへ入る。
既知boundary点を別sourceへ逆に用いた論理ではなく、全magを保つ実source写像そのものを構成した。

これでdouble_Aの全E義務を閉じ、新成分数を変更しない。
nil22 imageの未帰属はstable_line_conic、stable_three_linesの二源だけである。
他のE/Gの全支持、全分類、Y以外の未知親又は停止中の局所idealの精密支持は結論しない。
909主原稿/PDF・global台帳・Gitはこのreviewでは変更していない。
