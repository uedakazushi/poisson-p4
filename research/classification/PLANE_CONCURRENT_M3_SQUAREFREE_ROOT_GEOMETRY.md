# Concurrent三線・位数三SF：全root coversと有限first-connection入口

2026-11-08、G、producer plane_collision_review。純証明、独立review前。
新source probe、minor、GB、jet又はCASを行っていない。
Rは三相異なるconcurrent lines、f=a b(a−b)、第三base座標c、normal variables u,v。
πS rank二、actualD被約既約quintic、normal位数三、generic G3 squarefreeを仮定する。
全first connectionを有限algebraic constructionへ還元するが、実行可能profile listの全列挙とwhole higher G上界は未完了。
generic21/tripleの全完了は本稿の主張ではない。

## 1. 全root coversとc-degreeによる二分

U=(P1\{0,1,∞})×A1なのでroot monodromyは1,C2,C3,S3を全て許す。
任意のsplitting coverはbase P1のmarked cover curveからpullbackされる。
そのcurveのfunction field K上で、root vectorsをK[c]のprimitive polynomial pairsに取る。
pair determinantsとoverall coefficientは、全cでSFなのでK[c]のunits、すなわちK*である。
各root vectorのc-degreeの和はdeg_c G3≤2。
positive degreeが一つだけなら、他二constant independent vectorsとの二determinantsがconstantであるため、そのrowもconstantとなり矛盾。
従って全roots c-independent、又はdegree pattern(0,1,1)だけである。

connected C3/S3ではmonodromyが三rootsを巡回するため全c-degreeが等しく、後者は不可能。
split又はC2の後者では唯一c-independent rootはmonodromy-fixedでC(a/b)へdescendする。
global primitive factor ell(a,b;u,v)のbasedegreeが正なら、残quotientのbasedegreeは一以下でc²項を持てない。
ゆえにellは**constant normal root**、actual normal GL2でuへ直せる。
これで全SF源はc-independentと、constant uを持つc-dependentの二枝に尽くされる。
Galois coverのPic零やglobal UFDを仮定したのではなく、generic base field K上の一変数Gaussを使い、fixed factorだけをglobalへ戻した。

## 2. c-independent sourceはP1上のfinite-degree root maps

G3のcoefficientsは全binary quadratics in a,bとなる。
primitive root incidence curveのbidegreeは(e,3)、e≤2。
connected degree三root curveのnormalizationからnormal P1へのmap degreeはeである。
normal constant GL2はこのmapへのpostcomposition PGL2として実際に作用する。
overall base factorのdegreeは2−eで、そのdivisorはRに支持される。
e=1ならその全factorは三line componentsの一つ、e=2ならconstantだけ。
root-cover curveをrationalに変数変更することをambient PGL同一視とは呼ばず、
その上の全normal mapsを分類した後のpostcompositionだけをactual normal GL2として使う。

## 3. Constant normal rootを持つ全c-independent SF

split SFの全primitive degree flagsにはconstant normal rootがある。
C2でfixed root degree零の場合も含めて

`G3=u(Au²+2Buv+Cv²)`、A,B,Cはbinary quadratics、
`T=AC−B²`とCは非零、divisorsはR支持である。

全connectionは

```
p=dlogf−(1/2)dlogT,
q=dlogf−dlogC,
r=((B/2)dlogT−dB)/C,
β=[[p,0],[r,q]].
```

Cのdivisorはdegree二の全六monomials、Tはdegree四の全十五monomialsという有限集合。
nonzero scalar amplitudesはp,qへ寄与せず、rのcommon amplitudeはconstant normal diagonalで調整出来る。
B→B+tCはconstant normal unipotentであり、connectionをそのconjugateへ変えるだけ。
**全leading amplitudeを固定してしまう操作ではない。** 同一connection上の全parallel tensorsは後のsourceへ残す。

Cが二distinct linesの積なら、B²+T divisible Cの条件は二P1 pointsでのBの値を各有限なsquare-root choicesへ固定する。
binary quadratic Bは三係数で、残り一方向はCの倍数だけ。従ってB mod Cは有限個。
C=L²なら、ord_L T=0はBのvalueとfirst jetを有限に固定、ord一は不可能、ord≥二ではB=L(sM+tL)。
tはunipotentで除け、sはzero/nonzeroの二場合をnormal diagonalで扱う。
さらにglobal βのsimple-pole条件を全三componentsで検査し、未検査候補をactual connectionとは主張しない。
したがって全constant-root SFの**first connections**は、normal GL2を除いて有限algebraic constructionで得られる。
leading cubicsそのものが有限個との主張ではなく、各connectionの全parallel spaceを保つ。

## 4. C2、moving rational fixed root degree一の全有限候補

branch pairをa,bとし、cover curveをr²=a/bとする。
fixed root degree一、residual quadratic coefficient degree一なので、moving rootsのnormal mapはdegree一。
normal postcompositionでmoving mapsをr,−rとし、fixed rootはMöbius mapρ(a/b)となる。

`ρ(t)=(A t+B)/(C t+D)`。
ρ(r²)=rのcollision equationは

`C r³−A r²+D r−B=0`。

その全rootsは{0,∞,1,−1}に支持されなければならない。
degree三divisorをこの四pointsへ分ける全候補は高々20個で、各polynomialはρをscalar以外一意に決める。
AD−BC≠0を課し、全simple-pole/global conditionsを戻す。
ρ(r²)=−rの条件はr→−rで同じ有限supportに入る。
これで全C2 moving-fixed-root源を有限に覆い、arbitrary representative pointで置換していない。

## 5. C3 branch二lines：degree二involutionの全有限候補

base characterをinverseとline permutationsで割るとbranch pair a,b、third line a−b。
cover curveはP1_r、r³=a/b、deck σ(r)=ζr。
normal mapφのdegreeは一又は二。
degree一はnormal PGL2を除きφ=rで、primitive G3=a u³+bv³、overall hはRのline component。

degree二ならφの唯一involution τ=[[A,B],[C,−A]]を用いる。
φ(r)=φ(ζr)となる全rはbranch preimages

`D5={0,∞,1,ζ,ζ²}`

に入らなければならない。
τ(r)=ζrのquadraticは`ζC r²−A(1+ζ)r−B`。
全二degree divisors supported on D5は高々15個で、各polynomialはτをscalar以外一意に決める。
detτ≠0を課し、各τのdegree二quotient mapをnormal PGL2後の代表とする。
逆のcollision φ(r)=φ(ζ²r)はσで移した同じ有限supportなので、別coverを落としていない。
これが全degree二mapsの有限constructionである。

## 6. C3 branch全三lines：elliptic coverの全空性

cover curveはsmooth elliptic cubic C:w³=a b(a−b)、deck σはoriginをfixed pointに取ればmultiplication ζ。
normal mapのdegreeは≤2、elliptic curveからP1へのnonconstant mapはdegree二以上なので正確に二。
degree二mapはline bundle F∈Pic² Cのcomplete pencil、involutionはτ_s(x)=s−xである。
φ(x)=φ(σx)でσx≠xとなるpointは `(1+ζ)x=s` の唯一solution。
1+ζ=−ζ²はunitであり、U上collisionを避けるにはこのsolutionがσの三fixed pointsに入らなければならない。
従ってsはfinite fixed subgroupに入る。

H=O_C(1)はσ-invariant degree三。
H−Fはdegree一でunique point pを与え、上のconditionはpがramification三pointsのどれかという条件に等しい。
従って全φはこれら三pointsからのprojection pencilsである。
例えばp=(a=w=0)のpencil a,wのnormから、全leadingはbase permutationsを除き

`G3=a²u³+b(a−b)v³`。

そのdiagonal connectionではq=(1/3)dlogf（またはnormal swap後pがcanonical）。
v³ coefficientはfullhigherのconstant derivative C(f)の像に入らなければならない。
Q2=b(a−b)に対しC(f)=kQ2なら、mod b,mod(a−b)からC(b)=C(a−b)=0。
aもこの二formsのlinear combinationなのでC(a)=0、k=0となる。
coefficient非零に反する。従って**全elliptic C3 branchは空**。
この議論を他のrational coverへ無断で拡張しない。

## 7. Connected S3：全degree二mapを330個以下のdivisorsへ

Galois cover curveはP1_r、D3作用はr→ζr,r→1/r。
base coordinateはt=(r³+r^-3)/2で、branch valuesは−1,1,∞。
root fieldはinvolution r→1/rのfixed field、q=r+r^-1、t=(q³−3q)/2である。
normal root mapはψ(q)、degree一又は二。
degree一はnormal PGL2を除きq、primitive cubicは

`b(u³−3uv²)−2a v³`

であり、overall hはRの全line componentsを許す。

degree二ならψのinvolution τ=[[A,B],[C,−A]]を取る。
ψ(q(r))=ψ(q(ζr))の追加collisionはτ(q(r))=q(ζr)から得る。
denominatorを払うと

`Cζ r⁴+Aζ² r³−(C+B)r²+Aζr+Cζ²`。                    (2)

全rootsはbranch preimages

`D8={0,∞} union {sixth roots of 1}`

に支持される。q(r)=q(ζr)自身のcollisionもD8上にある。
degree四divisorを八pointsへ分ける候補は高々binomial(11,4)=330個。
(2)のcoefficientsは(A,B,C)へinjectiveなので、一divisorからτはscalar以外高々一つ。
detτ≠0を戻し、degree二quotient ψを取る。
全candidateの一つ一つをまだ算術列挙してはいないが、任意点を選ばず全coverを含む明示有限constructionである。
以後global β及びKが存在しない候補は捨ててよいが、未検査を空と呼ばない。

## 8. c-dependent SFは二つの全first profiles

§1によりfixed rootをu、G3=u(Au²+2Buv+Cv²)とする。
C,T=AC−B²はc-independent、R支持。
B=B_c(a,b)c+B0、B_cはnonzero linear。
polynomial Aのc² coefficientからC|B_c²。
従ってC=L²、B_cはL倍。mixed c coefficientからL|B0。
actual c-shear/scaleで全BをLcへ戻せる。
T=L²M2、M2はR支持のquadratic。
全connectionは

```
p=dlog(f/L)−(1/2)dlogM2,
q=dlog(f/L²),
r=(c/(2L))dlogM2−dc/L.
```

M2がLを含めばrはLに二位poleを持つので、R reducedのglobal条件に反する。
従ってM2は他二linesの平方、又は積だけ。
marked line permutationsで全sourceはこの**二first profiles**へactualに尽くされる。
M2のnonzero amplitudeはdlogへ寄与しないためleading全parallel spaceへ残す。

L=a、他lines b,a−bとして、raw P0の全literal fieldsを記す。
P_fは共通、{x,u}=uX、{x,v}=uZ+vY、normal K=K20u²。

| M2 | X | Y | Z | K20 |
|---|---|---|---|---:|
| b² | (a+b)∂c | (−2a+4b)∂c | (2a−3b)∂a−b∂b−c∂c | 2 |
| b(a−b) | (−a/2+b)∂c | (−2a+4b)∂c | (3a/2−3b)∂a−(3b/2)∂b | 0 |

K20はD2零からp(Z)、K11=K02=0である。
U上のrank二P0はflat horizontal connectionのbivectorなので、Kを戻したこれらはprojective Poisson。
flat接続だけからKのpolynomial性を推測せず、上表でconstant値を実際に解いた。

z=c u+a vとすれば全regular parallel G3 spacesは

* square型：span(b²u³,b z u²,z²u)、dimension三。
* nonsplit型：span(b(a−b)u³,z²u)、dimension二。

先頭coefficient一点を固定して全higher sourceを代表させない。

## 9. First-source countsと、まだ残る仕事

c-independent βはcを含まず、whole P0は∂c wedge W0(a,b,u,v)へ入る。
c-shearはこれをscalar以外変えないので、R choices五＋anchor amplitude一＋normal conjugacy高々三からB≤9。
c-dependent square型でもc-shift along the other root lineはconnectionを変えず、along Lはnormal unipotentへ吸収されB≤9。
nonsplit型はL以外の一c-shift directionが残り、安全なB≤10。
c-scaleはnormal diagonal、c-shift along Lはnormal unipotentと同じconnection変化なので重複して数えない。
fullhigher d≤10（B9）又はd≤9（B10）が得られればplane六/tensor scalar一を戻してG≤24が十分。

有限constructionの実profile列挙、wholeP0の全K/全350確認、全P1/P2とactualD3 realizabilityはまだ実行していない。
特に有限candidate constructionと全G完了を同一視しない。
known concurrent m三23componentの存在と両立し、全空性を主張しない。

## 10. Repeated02への共通canonical補題の独立適用

`PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md`、SHA
`4fd66aea9bda0adc00c5b4fe500fbac21b142a109da1aed709c9dd0d9ce34347`
を独立に手確認した。
P0 X,Yがc方向、K11=0、bad C=λ∂c非零ならcanonical divergence-zero liftは

`barX=(4X−Y)/5, barY=(4Y−X)/5, barZ=Z−(2K20/5)Ebase`。

P1のnormal divergence零によりbarC=C。
strict Jacobiのuv² coefficientは

`(2barK20 I+DbarZ)C+(DbarY)B−(DbarX)G=0`。

a,b projectionでは他項が零、Euler補正もC方向なので `(DZ C)_(a,b)=0`。
これはglobal stabilizerを仮定したcanonical化ではなく、全projective Poissonに一意なBV liftである。

repeated02 reduced determinant B=f/LでL=aとすると、全offdiag fieldはZ=3∇L×∇A2。
q=3dlogL−dlogB=dlog(L³/B)は三components activeなのでbadCは∂cに限られる。
∂cZのb componentは−3A_cc、従って全actual sourceに

`A_cc=0`, すなわち `A2=c ell(a,b)+Q2(a,b)`

が必要である。square determinantの全Aは既degree還元でc-independent。
これはrepeated02/triple中央wallへ使える新necessary conditionだが、その全higher G上界はまだ残る。
