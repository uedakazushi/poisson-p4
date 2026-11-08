# Smooth cubic scroll A25全source：二境界の同一定数と実軌道の独立監査

2026-11-08。PRIMARY E。別作者による**純証明PASS**。
対象は `CUBIC_SCROLL_A25_WHOLE_RESIDUE_CONSTANT_REDUCTION.md`
（7c286d08…）。新CAS、sourcepoint、operator、minor、profile又はjet無し。

受理する正確な結論は、A25 nodeの全actual源は空、cuspの全actual源は一つの
projective PGL軌道又は空、である。存在する一代表の既知親への実同定は別義務である。
そのactual同定を本receiptだけで認証したとはしない。

## 1. 全source、正しい座標、full-P2境界

既受理のscalar connection／quadric net／exceptional restriction一意性と、
`CUBIC_SCROLL_A25_BOUNDARY_SEPARATE_REVIEW.md`、座標・pole erratumを入力とする。
Y=Bl_S P4=P(F)→P2_(a:b:c)、E=2H−B、D′=H+2B。
全incidenceはa x3−b x2+c x0=a x4−b x3+c x1零。
a=x2一chartではx3=b−cx0、x4=b²−bcx0−cx1、
Q=x2x4−x3²=c(bx0−x1−cx0²)。したがってgeneric L=(c零)でord_L Q一。
旧符号誤記を再利用していない。

nodeの両root rows、cuspのdenominator rowだけをprimitiveとする。
cusp numeratorのdivisorial zero、D′のfull-P2 fibers、rootの重複、無限遠とRの全multiplicitiesを保持する。
c∤Ω、node22/31、cusp denominator double rootへの旧全境界還元を使う。
node31の既空性は再開せず、normalization-flat性を追加仮定しない。

## 2. 同じkとtotal exterior derivative

ordered rowsのLmat=[[-V1,U1],[-V2,U2]]又はindices(1,0)でdetLmat=Ω。
Sym²Lmatのflat coordinatesでQ′=Ω²Q、f=Y′/Z′、R=Q′/Z′²。
flat connectionとparallel vertical Poissonの全入口により、同じ一つの非零constant kで

    ωB=da∧db/Ω   (c一),
    ωf=(1/k)dlogR∧dlogf  (node),
    ωf=(1/k)dlogR∧df     (cusp)

となる。dはbaseも含むtotal derivativeである。
P2 chartのx=X′/Z′、f=Y′/Z′ではnodeのvertical Poisson coefficientはk(x−f²)f、
cuspではk(x−f²)なので、上式の相対符号とnormalizationを独立に確認した。
base chartや境界ごとにkを独立に再選択してはいけない。

a一chartではωB=db∧dc/(cΩ)でres_L ωB=−db/ΩL。
primitive denominatorとκ非零からZ′は各generic boundary fiberでunit。
ord_L R一、generic reduced Ωのc非零部分でord_Ω R二である。
したがってfiber residueはそれぞれ(1/k)dlogfL又はdfL、(2/k)dlogf0又はdf0。
fのtransverse derivativeはregularで、dΩ/Ωとの外積のpole係数はこの式そのものになる。

L及びΩのgeneric vertical divisorはactual E又はD′ではないので、全ωはそこでregular。
codimension二のfull-P2 fiberを除いたgeneric divisorで必要恒等式を得ており、
そのzero sourceを捨てていない。

## 3. Nodeの全正常化と二つの符号の矛盾

旧全22分類は両Aiを相異なるfixed pointのpure nilpotentとする。
actual PGL2と許された各row scaleでNlower／Nupperへ移す。
block shearの全r_i変換を手計算すると、残りはr1=(0,R0)、r2=(S0,0)。
両primitive性がR0S0非零を与える。
diag(d1,d2,k0)後のrow-scaleは
R0→k0R0d1/d2²、S0→k0S0d2/d1²。
d2一、d1³=S0/R0、k0=d1²/S0で両方一になり、全unitの逆と有限root選択だけで全sourceを覆う。
nodeの任意GL2 row mixingは使っていない。

canonical rowsは(b,c),(c,a)、Ω=ab−c²。
Lでa一、collapsed conic parameter z=bよりfL=−b²、ΩL=b。
従ってres_L ω=(-1+2/k)db/bで、k二が必要。
Ωでc一、a=1/b、f0=b、base residueはdb/b。
従ってres_Ω ω=(1+2/k)db/bで、k負二が必要。
同じ一つのkに矛盾し、node22も全空。
row交換でfを反転しても、両必要値が同時に符号を変えるだけなので矛盾は消えない。

## 4. Cuspの全denominatorとΩの補集合

A0=Nlower(n)、r0=(R0,S0)、nS0非零。
columns(e1,T0e1,T0²e1)のdet=n²S0なので、pを固定するactual共役でT0=N3に戻る。
全numeratorはold flag theoremによりlower triangularで、presentation row-addition後は

    (U0,V0)=(b,c),
    (U1,V1)=(λa+Pc, μb+Qc),
    Ω=λac−μb²−Qbc+Pc², μ(λ−μ)非零。

λ零をまだ捨てない。Ωにmultiplicity m≥2のcomponentがあれば、
baseωはorder m、fiberωはR=M^(2m)unit、f regularゆえ高々simple poleで、相殺出来ない。
したがってΩはreduced。
λ零ならreduced二linesで各f0がconstant、fiber residue零だがbase residueはgeneric非零。
よってλも非零で、Ωのdet=λ²μ/4非零からsmooth conicになる。
reducible Ωの段階で誤った「residueのpoleはL上だけ」を使っていない。

## 5. 全constant comparisonと残P,Qの実逆

L上a一ではfL=(λ−μ)/b、ΩL=−μb²だから
k=μ(λ−μ)。Ω上c一ではa=(μb²+Qb−P)/λ、f0=μb+Q、
base residue db/λ、fiber residue 2μ db/kなのでk=−2λμ。
全nonzero条件を確認した上で比較するとμ=3λ。
numerator scaleと対応するprojective tensor scaleでλ一、μ三へ戻る。

T0=N3、T1=[[1,0,0],[0,3,0],[P,Q,0]]とし、
H=I+uN3+vN3²でactual共役すると、off-diagonal entriesは
−2u、Q+3u、P+v−3u²−uQとなる。
presentation row-addition2uT0の後、Q→Q+5u。
u=−Q/5、v=−P+3u²+uQで全P,Qを零へ移し、
H逆=I−uN3+(u²−v)N3²、detH一であることも手で確認した。
P,Qのzeroを反転又は除外していない。

各全presentationがkerdfとσ、従ってπ_Eをprojective scaleまで一意に決め、
既global vanishing H0(Y,Λ²TY(−E))零が全ambient extensionを一意にする。
よって全cusp actual源は、(N3,diag(1,3,0))の一actual軌道又は空となる。
seedのみに基づくhigher推測でも、G次元上界からのE帰属でもない。
genuine一代表のactual T37同定を別receiptで受理した場合に限り、
この全正常化と閉PGL安定性からwhole A25 cuspをT37へ戻せる。
