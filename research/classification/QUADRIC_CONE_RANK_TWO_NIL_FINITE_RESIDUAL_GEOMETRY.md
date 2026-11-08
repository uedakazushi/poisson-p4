# Rank二quadric cone・nil：全既約D源に残る一whole族

2026-11-08、PRIMARY E。純手証明、新operator／profile／jet／CAS無し。
前の `NIL_BASE_REDUCTION` の全sourceを使う。
**唯一残す必要whole族はkappa=3/4,sigma=0、A任意で、
actualF=alpha*h^5−(5/2)c(q+2A*h²)²、alpha非零。**
その族の存在・全Y compatibility・E帰属は本稿ではまだ結論しない。

## 1. 通常のbase形と全Pfaffian

q=ac−b²、N=2b∂a+c∂b、w=1−2kappaとする。
kappa≠0,1では前noteにより全源の必要形は

```
B_ab=A+2kappa*q+b²,
B_ac=2bc, B_bc=c²,
T=sigma*Euler3+w*N,
sigma*A=0, T0=0.
```

全Y degree≤3、Y3=−q*Euler3を保持する。
column divergenceはZ_B=−2(1−kappa)N。
Faff=G+vK、G=bvec dot Y、K=bvec dot Tで、

`K=sigma*c*((1+2kappa)q+A)`。

全modular恒等式 `P gradF=−F divP` は前SS noteと同じ規約。
K非零では正確に `Y=gT−B#dg, g=G/K` となる。

## 2. sigma非零：全c-poleのnilpotent矛盾

sigma非零ならA=0。
kappa≠−1/2ではK=sigma(1+2kappa)cq非零。
Fにc因子がないと仮定すると、gのc=0でのsimple-pole residue g0(a,b)は非零。
c=0でbvec=(0,0,(1−2kappa)b²)なので、g0はY_c|c0の非零constant倍、
従ってa,bの多項式、degree≤3である。

B#dc=cNだから、Yのpole coefficientは

`g0*(T+N)−(B#dg0)|c0=0`。

そのb成分は正確に

`sigma*b*g0+w*b²*partial_a(g0)=0`。

w=0なら直ちにg0=0。w≠0ならC(b)[a]上で、非零constantと微分の和
`sigma+w*b*partial_a` は有限degree多項式にkernelを持たない。
最高a-degreeから順に全係数が零となるため、やはりg0=0、矛盾。
従って全sigma非零枝はc factorで空。

kappa=−1/2ではK=0、Faff=G。
modular恒等式のc成分はN(G)=−sigma*Gを与える。
Nは多項式環上locally nilpotentなのでsigma非零に対応する固有多項式は零。
従ってここでもG=0で対象rank四でない。

## 3. sigma零：全N不変polynomialのweight比較

sigma=0ではK=0。
modular恒等式のc成分からN(G)=0。
標準nilpotent作用のpolynomial kernelは正確にC[c,q]である。
例えばcを反転してN=c∂bと直し、c-adic poleをq=ac−b²で比較すれば
intersection C[c,c^-1,q]∩C[a,b,c]=C[c,q]を得る。
各c-adic最低係数でq modc=−b²の独立性を用い、c零のsourceを捨てていない。

B#dc=cN、B#dq=(wq−A)Nなので残る式は

`[c∂c+(wq−A)∂q]G=(1+2w)G`。

w非零ならQ=q−A/wとおき、全polynomial解のmonomial c^p Q^jは

`p+w*j=1+2w, p+2j≤5, p,j≥0`

を満たす。c因子を避けるにはp=0の項が必要だから、j≤2よりw=−1/2又は−1。

* w=−1/2、kappa=3/4：全解はalpha+beta*c*Q²。
  固定Y3からbeta=−5/2、Q=q+2A。
  alpha零はfactor、alpha非零のwhole族を残す。
* w=−1、kappa=1：全解はalpha*Q+beta*c*Q²でQ因子。

w=0ではc^p*g_p(q)ごとに `−A*g_p'=(1−p)g_p`。
有限degree多項式ではp=1だけが残り、全Gはcで割れる。A零も同じ結論である。

## 4. kappa=1の全補集合

前のbase/T必要還元により、ここでは既に

`B=Jac(c(q+A)), T=sigma*Euler3−N, sigma*A=0`。

これは§1と同じB,Tの形なので、sigma非零には§2、sigma零には§3をそのまま適用できる。
すべてproperfactor又はG零で空。kappa=1に別のA,d,D sourceを隠していない。

## 5. kappa=0の全(r,s,A,T0)を残す

前base noteのt=d=D=0を入れると

```
B_ab=A−rb+sa+b²,
B_ac=c(2b−r), B_bc=c(c−s),
T1=sigma*Euler3+N.
```

Lie_T B=0と[T,Z_B]=0を全係数で比較するとT0は

`T0=(tau,−sigma*r/2−s,−sigma*s)`

で、さらに `sigma*s=0, s(r+tau)=0` が必要である。
以下はこの全zero／unit分岐を扱う。

### 5.1 sigma非零

s=0、b'=b−r/2でBは§1のkappa0、A'=A−r²/4へ戻る。
Tはsigma*Euler3+N+delta∂a、delta=r+tau。
actual a'=a+delta/sigmaでlast constantを消す。
Bはaに依存せず変わらない。Lie_T Bのconstant係数がsigma*A'=0を与える。
したがってA'=0で§2のpole矛盾が適用でき、c因子で空。
すべてactual affine base変換で全P4へlinearに延び、全Yもpushする。

### 5.2 sigma=0、s非零

tau=−r。actual b'=b−r/2,c'=c−sとおく。
T=N'=2b'∂a+c'∂b'、q'=ac'−b'²、A'=A−r²/4。
全Bは正確に

`B=B2(kappa0)+s*Jac(q')+A'*∂a wedge ∂b'`。

従ってZ_B=−2N'、K=bvec dot T=0、B#dc'=(c'+s)N'、B#dq'=(q'−A')N'。
modular恒等式はN'(G)=0を与え、全GはC[c',q']内。
残りは

`[(c'+s)∂c'+(q'−A')∂q']G=3G`。

U=c'+s（元のc）、Q=q'−A'とすると全解はtotal degree三のU,Q homogeneous polynomial。
元の空間次数が高々五なのでQ³項は無い（その最高次数六は他項で相殺出来ない）。
従って全項にU因子があり、元c因子で空。

### 5.3 sigma=s=0

b'=b−r/2によりB=B2(kappa0)+A'∂a wedge ∂b'、T=N'+delta∂a、delta=r+tau。
delta=0なら§5.2のs=0版でc因子。
delta非零ならK=delta*c²。
c因子を避けるGはG|c0非零なので、g=G/Kはc=0にorder二のpoleを持つ。
そのleading coefficient g0(a,b')はG|c0/delta、多項式degree≤5である。

Y=gT−B#dgのorder二pole全体は

`g0*(T+2N')−(B#dg0)|c0=0`。

b'成分から(b'²+A')∂a g0=0、従ってg0はb'だけの多項式。
a成分は

`(b'²+A')*g0'=(6b'+delta)*g0`。

非零polynomial g0の最高次数をnとすれば、最高係数比較がn=6を強制する。
しかしn≤5なので矛盾。A'の零／非零又はquadraticの重根を一切反転していない。
従ってc因子で空。

これでkappa0の全unit／zero補集合を尽くした。

## 6. 残る義務

constant-C nilpotentの全既約quintic源は、§3の唯一のwhole族へ入る必要がある。
残るY compatibility、alpha非零での実存在、全normalizerとactual既知parentを次に扱う。
本稿から新componentや全E完了を主張しない。C零枝の既空性、SSの別有限入口と混同しない。
