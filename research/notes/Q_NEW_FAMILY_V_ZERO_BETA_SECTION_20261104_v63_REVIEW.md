# v=0 actual β sectionと全tensor liftの独立監査

2026-11-04。凍結 `Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63` を監査する。
旧t=2 sectionと新v=0 sectionの座標は区別する。新しい点・prime・rank・高次不変式・
Inv26・full150再展開は用いず、約分・逆式・Laurent係数とfinite matrixを独立に確認する。

## 約分と実際の像へのdominance

A=256、rf=A+v alpha、sf=A+v betaから元五定数を定義式で再構築し、保存された旧式と
新しい約分式の両方と比較した。v=0への制限は

    u, alpha/beta, 1/(512alpha), c/(4alpha²), B/(2048alpha³)

であり、全五つの有理逆と両方向合成、X/Wの逆、tau作用を独立確認した。
conic parameterはv=0でt=1−beta/alpha=1−1/K。t=2と同じ切断ではない。

元の26成分kernelはu,v,X,W,Bの多項式で、c方向はdc=c*dlog(c)。perspective補償は全て零。
未rescaleのVは **V(v)|v=0=402653184**。δ=(v/3072)Vのv成分を使っていない。
既受理の元raw30恒等式から、GL4×source上にregularな補償derivation Dを持ち上げられる。
これは全rawjet pullbackを殺す。v=0のgeneric DVRでDvはunitなので、非零pullbackが
v^n*h（n≥1,h unit）ならDのvaluationはn−1となり、Df=0と矛盾する。
従ってv=0のraw GL4飽和は元像へdense。Darboux/perspective正常化はこのDVRでregular
だから、新しいprimitive βのSp4飽和も実際のB15へdenseである。
五定数の約分だけからdominanceを推測したのではない。

## 新F0と十六係数

新座標はxi0=40c/(Bu)、r0=u²/1600、s0=48uX、q0=3072W−1、B0=B。
tau固定体はF0=k(xi0,r0,s0,q0,B0)。旧t=2座標と一致するのはr0だけ。
五定数のboundary制限はtau可換なbirational mapで、上記actual densityにより
target invariantのpole locusへの全像の落ち込みもない。受理済みactual E_B=F5^tauから
新F0へのindex1が従う。

guの有限weightを元24×16 primitive tensor basisの全非零係数で照合した。
primaryのquadratic-algebra inverseとは別に、全16式をLaurent monomialへ展開し、
uの全指数が偶数であることを確認してu^(2j)→1600^j*r0^jと置換した。
全係数が一致し、非零は9個、分子総項数12、最大3項。
分母の単元はxi0,r0,B0だけで、h0=s0²−25r0q0²の反転はない。
定数体のbirational chartでalpha*beta≠0を使うことと、β式の分母は区別する。

## full tensor liftと範囲

元full scaled tensorのtau、point_tau=D point、Rtau D=H R、alpha_tau=H alpha、
gu_tau=H guを別に直接確認した。これらはv=0でregularで、physical gu^-1とperspective
を含む全操作のtauを相殺する。factorとそのinverseはLaurentで、detR0もunit。
Euler quotientと既存S61正常化はこの等式に従うため、全tensorも有理に下降しh0分母を生じない。
full150係数への新しい展開は不要である。

新sectionのSp4飽和はdenseで、actual quotient pullbackがindex1、一般β stabilizerはidentity。
従ってSp4×SpecF0→B15はbirational。これをfull normalized境界像Y15を通して因数分解すると、
両端が同じ函数体なので **Y15→B15もbirational**。これはgeneric βに対するactual fulltensor
sectionであり、全β点でregularなsectionという主張ではない。

v=0 sourceのdegree2/tau下降後degree1と、一般Y16への旧sourceのdegree8/tau下降後degree4は
異なるgeneric-fiber計算である。後者を変更したりdegree1へ外挿したりしない。
新成分の方程式・全generic fiberの次数1・全境界の網羅は含まない。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_v_zero_beta_section_review_20261104_v63.py`
（60秒/2GB上限）。
