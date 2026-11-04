# Reciprocal四次商・全二次cover・種数上界17の独立レビュー

2026-11-04。原著RECIPROCAL_QUOTIENTとGENUS_BOUNDの二packetをまとめて監査した。
source/β点、specialization、GB、octic因数分解は使わない。

独立したFraction/Laurent係数演算で、元の凍結Q1,Q2が正確に `RT*L1,RT*L2` であること、
reciprocal involutionが両式を保つこと、r,t,zが固定・hが反変であることを検証した。
普遍的なquartic identity、逆座標の分子identityも全係数で確認した。
実パラメータc,d,f,u,v,K1,B0は、原著のcleared bar式とは別に**直接Laurent式で代入**し、
全高次相殺と1592整数係数を照合した。plane式は `16*x²*y^10*Psi`、全次数4、
r4係数16x4y8、quadratic coverのqは全次数3である。
4×4行列でも同時負号・共役・diagonal normalizationを別に検査し、reciprocalの向きと
off-diagonal符号が一致した。

## 全generic algebraとしての表示

generic baseは x,y,x²−y²,m,c,f,u,v を単元とする凍結ordered fieldである。
quartic base `A0=K(m)[r]/Psi` はleading coefficient x²/y²が単元なのでfree4。
q=r²+4fuvがその中で単元であることは、既約性を仮定せずに証明できる。
もしPsi=qの共通幾何根があれば、rはk1を除く底上代数的だが、
t=(K1−xr)/yはk1に一次で超越的。W=0の定数係数からB0+2ar=0、
従ってB0²+16a²fuv=0が必要となる。しかしそのm4係数は独立計算でも4であり矛盾。
有限一変数algebraで共通根がないことはqの単元性を意味し、非被約でも有効である。
具体的なq逆元は4×4 multiplication matrixのadjugateで有限に構成できる。

h²=qを加えたalgebraはfree2、h^-1=h/qである。
R=(h−r)/(2v)、T=(z/h−t)/(2c)に対して、
`1/R=(h+r)/(2fu)`、`1/T=(z/h+t)/(2uv)` も得る。
両方向のgenerator合成がidentityなので、元Laurent idealと新idealは**全scheme algebra**として
同型。零点集合や選択したreduced factorだけを比較したものではない。
char0で2hが単元だからgeneric二次coverはétaleで、h→−h固定algebraは正確にA0。
special m・projective boundaryでの分岐やcomponentの選択は別問題である。

## 各幾何的componentの上界

係数の代数閉体上、mに優越する任意のcomponentを取る。対応するplane componentのdegreeは
d≤4、normalizationのgenusは `(d−1)(d−2)/2≤3`。
qは(r,m)でdegree≤3のpolynomialだから、normalizationへのline-at-infinity pullbackの
3倍以下しかpoleを持たず、pole divisor degree≤3d≤12。
zero divisorも同degreeなので、odd valuationによるquadratic ramification degree≤6d≤24。
connected二次coverではRHにより `g≤2*3−1+24/2=17`。
splitなら各componentはbaseと同じfieldでgenus≤3。

plane式がreducibleまたはnonreducedでも、各integral componentのnormalizationに同じ議論を
適用できる。四次の滑らかさ、既約性、正確なgenus、selected component degree8、
Weyl quotient前後のgenus一致は仮定も結論もしていない。
対象はordered eigenbaseと既存relative-sign quotientの後の曲線に限る。
