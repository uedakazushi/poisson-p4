# P*の真T19帰属、全5H源、六候補への全支持還元：独立査読

2026-11-08。integrator `/root`。producer `/root/generic_exhaustion` の固定実曲線と全h5合成を読解し、元true completed T19から同じ実曲線の全150/350/126を独立に再生した。新しい点、frame parameter、period、ambient rank又は局所変形jetを探索していない。旧限定criterion FALSEのreport/caller/jobを変更しない。

## 1. 一つの固定actual曲線とその逆

source(a,b,c,d,H)、target(h,x,y,z,w)。r=-e³,v=2e²-e³,t0=-e²,t1=-2e²。元T19のcoupled correctionを全て保持したcompleted homogeneous tensorを読む。eunitでperiodopen -r(r-v)=-2e5、actuallinear mapは

```
h=-2e³(b+c), x=e⁴(a-2b-4c), y=-e(2b+c+d),
z=e²(a+2b+2c), w=-a+2b+c+d+H.
```

det=4e10、逆は

```
a=h/e³+z/e²,
b=x/(2e⁴)-3h/(2e³)-z/(2e²),
c=h/e³-x/(2e⁴)+z/(2e²),
d=-y/e+2h/e³-x/(2e⁴)+z/(2e²),
H=w+y/e+z/e²+h/e³.
```

独立再生で両合成を確認した。tensor scale一、source/targetの全Euler gaugeを保持し、raw liftの負冪を消したと主張しない。canonical projection P↦P-E∧div(P)/5の全150coefficientsがe多項式となる。全350元Jacobiとstrict projective Jacobiは零。signedEuler-Pfの全126係数はdet共変則

```
D=4e10[-r(r-v)H(cd-b²)(ab+dH)](inverse)
```

に一致し、e0でh5となる。全150中心は以下のP*のcanonical sectionそのもの：h-row零、順(xy,xz,xw,yz,yw,zw)で

```
(hx,-x²,xz-h²,xy-h²,-yz+hw-hy,-xw+xy+hy+hz).
```

中心非零なので曲線はprojective parameter spaceへ正則に延長する。generic eは名前付き真T19のactual像、T19は閉PGL不変だからP*∈T19。

## 2. 旧FALSEと固定符号の監査

保存universal report32dd5437の全canonical150へ(Rpar,V)=(-1,2)を入れ、S=diag(1,1,-1,-1,1)でpushしたものが、上の独立再生の全150と完全一致した。S^-1=S,detS1。
旧保存中心は-V Rtilde/2-Rpar B、Rtildeのxw=xz+h²,zw=-xw-hy。全六置換でS_*Rtilde=-R,S_*B=B。従って(V/2)R-RparB=R+B=P*。
旧自動criterionは符号をpostcomposeせず中心をaR+bBに直接入れるものだったためFALSEである。新しいactual帰属はそのFALSEを成功へ改竄したものではない。独立再生は保存したuniversal全150との同一性まで検証しており、別curveを事後探索したものでもない。

独立exact report `CLASSIFICATION_UNIMODULAR_POLE_T19_FULL_SAME_CURVE_ROOT_REVIEW_20261108.json` は2.827372秒、最大RSS54360KiB、全チェックPASS。実行前root契約とinline-source callerを併存し、full150/350/126と正逆、rawEuler場を保存した。

## 3. 全h5への合成

既別作者受理 `UNIMODULAR_H5_WHOLE_ONE_ORBIT_SEPARATE_REVIEW.md` は全C6/L24/Q60とunit constant Pfの全解をknown閉和＋P*一軌道に還元した。六saturated highest closures、全lower fiber、gamma共通factor/zero、non-split/nil、nonsaturated thetaのdegree1/2 factor properquartic-D入口を別々に戻した結果であり、最高部のgenericclosureから全lowerを推測しない。
唯一alignedzero源はn1n2unitで、元全六からactualdiag/tensorとその全逆によりP*へ。全translations/highestGL4を戻してもT19はPGL不変。Pf0補集合はこのrank四h5源の仮定外で、rank二既定理へ分ける。
従って全A=c=0、Pfconstant unitのquadratic源がknown閉和に入る。actual D=5HではHはPoisson hyperplane、一意h-Casimir liftのsymplectic体積からdiv0となるので、これが全rank四5H locusの点別被覆である。
本runで全六highest分類や以前のnormalizer算術を再実行したとは称さない。新しい独立再生は§1–2の最後の同じcurveだけ。全h5合成の仮定・補集合は別途読解した。

## 4. 全被約支持の上界148

既受理global149 theoremは、全rankのP_redをknown142と七既約閉候補で覆う。§3により、そのうちunimodular P*の軌道閉包は既知T19内なので有限和から除ける。残る六つは

1. Ga9 linear-w一軌道閉包（dimension≤23）、
2. 二simpleH＋nodal cubic plane全A6源閉包（≤23）、
3. 二simpleH＋cuspidal cubic/flex配置 plane全A6源閉包（≤23）、
4. 二simpleH＋cuspidal cubic/tangent-secant配置 plane全A6源閉包（≤23）、
5. constant modular highest-double強制軌道閉包（≤23）、
6. constant modular highest-triple強制軌道閉包（≤22）。

有限閉和の極大既約部分集合は一項に等しいから、全成分数nは142≤n≤148。未採録成分はこの六項のどれか、dimension≤23。各六候補が互いに異なる／known親の外／全成分という主張はこの帰結に不要であり、ここでは証明しない。
新成分加算零。exactcount、残六の親・包含・成分性、主分類完了は未決。全非被約schemeやP*の局所smoothness／枝数も非主張。旧149のimmutable証明は弱いが正しい記録として保持する。
