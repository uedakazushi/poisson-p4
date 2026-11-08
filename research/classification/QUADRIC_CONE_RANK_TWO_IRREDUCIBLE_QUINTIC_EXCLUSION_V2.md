# Rank二quadric-cone全排除 V2：全Fu-only入口と追加例外を復元

2026-11-08、PRIMARY E、producer純証明。独立review前。
旧whole稿は `ENTRANCE_ERRATUM` により非零Cの入口欠落を明示し、bytesを保持した。
本V2は `FULL_BASE_REDUCTION_V2` の**元全33係数からの正常化**を使う。
新CAS／operator／point profile／jetは一度も実行していない。

**定理候補。** 被約既約actual quintic Dを持つrank四Poisson構造には、
一般restriction rankが二のintegral quadric-cone非正規支持は無い。
全Fu-only、全Y,T0、全kappa/sigma、全normalizer unit／zero補集合を含む。
採用は独立全chain監査後にrootだけが行う。

## 1. 不変な全source入口

元PENCIL全33／全350／全126と正しいOre帰結は旧独立全算術receiptをそのまま使う。
chart h=1でP=B+(Y+vT) wedge ∂v、Y degree≤3、Y3=−q*Euler3、q=ac−b²。
actual Faff=G+vK、G=bvec dot Y、K=bvec dot T。
全modular恒等式は `P gradF=−F divP`。
K非零なら `Y=gT−B#dg, g=G/K`。

constant C零は `ZERO_C_ACTUAL_FACTOR` と新独立C0 receiptにより、
Fv零H因子／非零の両orbit二次因子で全空。
この枝には誤って落とした非零CのFu-only追加base項は無い。

## 2. SS：全追加項を戻した後のcanonical枝

`FULL_BASE_REDUCTION_V2` §1.1はkappa≠−1の元全sourceをactual affine変換で

```
B_ab=−ab, B_ac=d−w*q−2b², B_bc=−bc,
T=rho*Euler3−wH, w=1+2kappa, rho*d=0
```

へ戻した。全Yはpushし、degreeとleading termを保つ。
したがって旧 `SS_FINITE_RESIDUAL_GEOMETRY` §§1–4のpole／monomial計算は、
ここで初めて正当な全source lemmaとして適用出来る。

具体的にrho(2kappa−1)非零ではK=rho(2kappa−1)bq。
b-pole coefficientはY_b|b0のconstant倍なのでpolynomial degree≤3。
その全tangent derivativeを含む式から
`w(i+j−2)=2, rho=w(1−j)+1`を得て、残るのはkappa=−3/2,rho=±1だけ。
rho零ではH(G)=0、Q=q−d/wにおけるweight比較から
kappa=−3/4のG=alpha+(5/2)b(q+2d)²だけ（alpha非零）を残す。
w零とkappa=1/2の全rhoはproper factor又はG零。
すべてpolynomial pole residueと全nonnegative exponentを使い、normal部分だけの不十分なpole比較ではない。

残二typeは旧whole稿 §§2–3の独立した代数計算で閉じる。

* kappa=−3/4：bvecはproper ideal(b,q+2d)に入り、G mod ideal=alphaなのでalpha=0、矛盾。
* kappa=−3/2：q-pole coefficient eta=(gradq dot Y)|q0/(4rho)はgraded ring C[a,b,c]/(q)のpolynomial。
  B|q0=Euler3 wedge H、pole cancellationはEuler3(eta)=−4etaを与えるのでeta=0。
  従ってFにq因子がある。

## 3. SS unimodular枝の全追加例外

kappa=−1はV2 base §1.2を使い、全f2を保持した。
xi=nu=0へ戻る全canonical cubic枝では旧finite pole論により、
rho=−1/2のf=−bq+c、rho=+1/2のf=−bq+aだけが残る。
これらは旧whole稿 §4の**有効なconditional field lemma**で排除される。
例えばf=c(1−ab)+b³でK=3rho*f、eta=G/(3rho) mod f非零を仮定すると
B#deta=eta*T、field C(a,b)上でeta²(1−ab)=constant非零。
prime1−abのodd valuationに矛盾する。反射caseは1−bcを用いる。
その他の全canonical cubic parameterは旧pole又はQ因子で空。

省略してはいけなかった二例外は以下で新たに閉じる。

### 3.1 sigma=2、xi非零

全actual affine正常化後

`f=lambda*a²+b³−abc, B=Jac(f), T=(3a,2b,c), lambda≠0`。

Tf=6f、K=6f。
Fにf因子がなければeta=G/6 mod fは非零のregular polynomial classであり、
pole cancellationがB#deta=eta*T on f0を与える。
fはprimitive c-linearでirreducible、fieldはC(a,b)。
B_ab=f_c=−abから

`a*partial_a eta=2eta, b*partial_b eta=−3eta`、

従ってeta=C*a²/b³、Cは非零constant。
しかしこれはaffine surface f0でregularでない。
c=1,b=tとして、quadratic equationのformal branch

`a=t*(1+sqrt(1−4lambda*t))/(2lambda)`

を取る。aはvaluation一、bもvaluation一なのでa²/b³のvaluationは−一。
全regular polynomial classはこのC[[t]] branchで非負valuationを持つので矛盾。
これは高次Poisson jet計算でなく、一つの明示algebraic branchのDVR valuationである。

### 3.2 sigma=−2、nu非零

f=lambda*c²+b³−abc、T=(−a,−2b,−3c)、K=−6f。
field C(b,c)上の同じpole equationからeta=C*c²/b³。
a=1,b=t,c=t*(1+sqrt(1−4lambda*t))/(2lambda)のbranchでvaluation−一となり矛盾。
両lambda零は既canonical枝へ戻っている。

これでSSの全追加Fu-only方向も含む全sourceが空となる。

## 4. nil：modular fieldをcenterした一般枝

V2 base §2.1（m=1−kappa非零、xi=0）の全sourceは

`B=B2(kappa)+ell*Jac(q)+h*c*(∂a wedge ∂b)+A*(∂a wedge ∂b)`、

T=sigma*Euler3+wN+tau∂a、w=1−2kappa、
sigma*ell=0、2kappa*tau−sigma*h=0、ell*tau−2sigma*A=0へ入る。
N=2b∂a+c∂b。

sigma非零はactual a-translationで全追加項を消し、旧canonical B2,Tへ戻る。
旧 `NIL_FINITE_RESIDUAL_GEOMETRY` §2のpolynomial residue nilpotenceがFにc因子を与える。
kappa=−1/2のK零ではN(G)=−sigma*G、locally nilpotent性からG零。

### 4.1 sigma零、kappa非零

actual a-translationでhを消し、T=wN、B=B2+ell*Jac(q)+A∂a wedge ∂b。
全K=0、GはC[c,q]に入り、

`[(c+ell)∂c+(wq−A)∂q]G=(1+2w)G`。

U=c+ell、w非零ならQ=q−A/wとおく。
全monomial U^p Q^jはp+w*j=1+2w、p+2j≤5。
U因子を避けるとkappa=3/4、又はkappa=1だけ。
kappa=1は§5で同形を正当化して扱う。
kappa=3/4ではG=alpha−(5/2)U(q+2A)²。
全bvecはproper ideal(U,q+2A)に入るので、G mod ideal=alphaからalpha=0。
idealの商はC[a,b]/(−ell*a−b²+2A)で全ell,Aに非零。零補集合を捨てていない。
w零もpolynomial coefficient比較で全GにU因子。

### 4.2 sigma=kappa=0

ell*tau=0を残し、B=B2(0)+ell*Jac(q)+h*c∂a wedge ∂b+A∂a wedge ∂b。
tau=0ならK=0、N(G)=0、

`[(c+ell)∂c+(q−A−h*c)∂q]G=3G`。

U=c+ell、V=q−A+h*ellでoperatorはU∂U+(V−hU)∂V。
全polynomial固有解はtotal U,V-degree三で、nilpotent部分hU∂Vが零。
h非零ならGはU³のconstant倍。
h零ならdegree≤5によりV³項が無く、やはり全項にU因子。

tau非零ではell=0、K=tau*c²。
c因子を避けるGのorder二pole coefficient eta=G|c0/tauはpolynomial degree≤5。
pole equationのb成分からetaはaに依らず、a成分は

`(b²+A)eta'=(6b+tau)eta`。

非零polynomialにはdegree六が必要で、degree≤5に矛盾する。
h*c項はc=0で消えるので、この議論は全hを含む。

## 5. nil unimodular kappa=1の全Fu-only

V2 base §2.2によりsigma非零はB=Jac(cq)、T=sigma*Euler3−Nへactualに戻り、§4のsigma非零pole論で空。
sigma零は

`f=(c+ell)q+A*c, B=Jac(f), T=−N`。

これは§4.1と同じB2(1)+ell*Jac(q)+A∂a wedge ∂bの全形である。
U=c+ell,Q=q+A、weight方程式p−j=−1とdegree≤5から
G=alpha*Q+beta*U*Q²、従ってproper Q factor。
全ell,Aとsecond-order sourceを含む。

## 6. nil kappa=−3/2、xi非零の新例外

V2 base §2.3によりsigma=0、T=−(4/5)Z_B、K=0、Z_B(G)=0が全sourceで必要。
actual affine変換後Z_B=(−10b,−5c,40xi)。
homogenizeしてhを追加し、非零constant rescalingで標準Jordan4

`D=h∂x1+x1∂x2+x2∂x3`

とする。degree五kernelの完全六基底は

```
h^5, h³I, h²J, h*I², I*J, h*Delta,
I=2h*x2−x1²,
J=3h²*x3−3h*x1*x2+x1³,
Delta=(J²+I³)/h².
```

Deltaはdegree四polynomial（分子がh²で割れることは直接展開で分かる）。
完全性はSym^5(Sym^3 C²)のhighest weights15,11,9,7,5,3から従う。
対応irreducible summand dimensions16+12+10+8+6+4=56で全dimensionに一致し、
六基底のweightsは相異なるので独立。これは既受理cubic-cone B_nil proofと同じ有限表現論である。

h=0で残るのはI*Jの−x1⁵だけ。
したがって元affine degree≤5 invariantのdegree五部はc⁵のconstant倍に限る。
しかし固定Y3と元B2から、全actualGのdegree五部は2c(ac−b²)²でa²c³係数が非零。
actual affine変換と非零diagonal rescalingでこの非零係数は消えない。
矛盾するのでこの新whole例外も空。

## 7. 全補集合と受理範囲

C零／非零SS／非零nil、全kappa/sigma、全Fu-onlyをV2 base normalizerへ戻した。
元sourceを必要supersourceへ広げた段階は排除に安全で、全actual inverseも示した。
rank四Pfaffianがproper二次／三次因子で割れる全枝とF零枝を合わせ、定理候補を得る。
chart因子は元actual degree五にdegreeを合わせてhomogenizeし、Euler gaugeでも変わらない。

旧17-file稿を修正なしで保存したこと、誤入口の発見とC零だけの独立受理、
今回新たに処理したsigma±2／kappa−3/2,xi非零を明確に区別する。
新componentの増分は零。全global登録はV2独立review後のroot判断に委ねる。
