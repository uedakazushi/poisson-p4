# constant DOUBLE は全 T36/T38 の外：別作者の純幾何照合

2026-11-08。担当 `/root/plane_collision_review`。
root の候補証明を、受理済み actual tensor と true completed T36/T38 の全式から独立に手展開した。
新しい point/operator/minor、局所変形 jet、GB 又は CAS を実行していない。
ここでの二次最低項は固定 actual polynomial の手展開であり、高次補正を探索する操作ではない。
結論は **DOUBLE の一つの全 actual orbit が whole closed T36/T38 に入らない**こと。
新成分性、Ga9/T37 帰属、全 prior23 の監査、又は全E完成は主張しない。

## 1. 対象と全 intrinsic zero support

座標(h,t,x,y,r)、h-row零とし、全十括弧のうち六つは

```
Ptx=h²+xh−2xy, Pty=4y²−yh, Ptr=(h−6y)r+xh,
Pxy=0, Pxr=−x², Pyr=h²+2xy.
```

これは全 highest-double 源の残余を一 actual orbit へ送る受理済み finite reduction の代表である。
h非零 chartでは Pxr=−x² と Pyr=1+2xy により projective rank零はない。
ここでは h が Casimir なので、h非零上の raw spatial brackets が intrinsic brackets と一致する。

h零上で Euler を掛けた四 cubic は、順序(txy,txr,tyr,xyr)で

```
−6xy², x(−tx+4yr), 2y(tx+5yr), 3x²y.
```

最後から xy=0。x零なら y²r=0、y零なら tx²=0 なので全零台は厳密に

```
L1={h=x=y=0}, L2={h=x=r=0}, L3={h=y=t=0}.
```

三本は L2—L1—L3 の chain で、support dimension は一。
L2\L1 では y非零で第一 cubic の x一次係数が −6y²、
L3\L1 では x非零で第四 cubic の y一次係数が 3x²。
従ってそこでは intrinsic section の first jet が非零である。
零点で Euler wedge は tangent quotient への fiberwise injection を与えるため、
cubic のこの一次項は projective bivector の一次項の非消滅を意味する。
以下の chart式で L1 の全点が order二以上であることも確認する。

## 2. 全 order二零点の最低項

r=1 chart、L1上の点 t=t0 で τ=t−t0、normal base を(h,x,y)とする。
projective affine bracket は Pij+xi*Pjr−xj*Pir である。
全六 entry の手展開は

```
πhτ=h²−6hy+xh²,
πτx=h²+4xy−(t0+τ)x²−x²h,
πτy=10y²−2hy+(t0+τ)(h²+2xy)−xyh,
πhx=−hx², πhy=h³+2hxy, πxy=xh²+3x²y.
```

よって全0/1jetは零で、二次最低項は Q=∂τ∧Vt0、

```
Vh=−h²+6hy,
Vx=h²+4xy−t0*x²,
Vy=10y²−2hy+t0*h²+2t0*xy.
```

残る一点 [t:r]=[1:0] は t=1 chart を使う。
二次最低項は Q∞=∂r∧V∞、

```
V∞=(0,x²,−h²−2xy).
```

これは h,r,x,y の全 projective chart correction により直接得る。
従って全 L1 を覆い、有限 t0 のみに制限していない。
全 Vt0 は Euler E3=(h,x,y) の linear multiple ではない。
実際 h成分から multiplier は −h+6y であるが、x成分には h² があり一致しない。
V∞ も h成分零より非零 radial multiple ではない。

## 3. 全点に共通する irreducible invariant cubic

```
C=x²y+xh²+h³.
```

C は yについて primitive linear で、gcd(x²,h²(x+h))=1 より絶対既約。
全 t0 で直接微分すると

```
Vt0(C)=(−2h+18y)C, V∞(C)=0.
```

例えば t0の全項は
−t0*x²*(2xy+h²)+t0*(h²+2xy)*x²=0 と相殺し、
残りは −2hx²y−2xh³−2h⁴+18x²y²+18xh²y+18h³y である。

四次元 tangent space 上の cubic Euler-horizontal one-form
αQ=i_(E4 wedge Q)vol4 は、base三変数上の i_(E3 wedge V)vol3 の pullback に scalar を掛けたもの。
上の非radial性より各 αQ は非零。
しかも αQ は irreducible cubic cone {C=0} を保つ。
Cが全 cubic coefficients を割るなら αQ=C*η、ηは constant one-form となる。
Euler horizontality は C*η(E4)=0 を与えるので η=0、矛盾。
従ってこの invariant hypersurface の一般点で foliation は非零である。

## 4. T36 の proper closed 条件との矛盾

この節は `GA9_LINEAR_W_WHOLE_T36_QUADRATIC_PENCIL_EXCLUSION.md` §§1–3 の
全source/properness 論証を独立に読み直したもの。
true completed T36 の順序(A,B,C,u,v)、k=p0/5≠0 とする。
q+ = kBu+t2*v²、q− = kCv−t3*u² は独立 quadrics である。
A=1 の vertex の quadratic jet は全 t0,t1 を保ったまま

```
QBC=5kBC+2t2t3*uv/k, QBu=3kBu+t2v², QBv=kBv,
QCu=−kCu, QCv=−3kCv+t3u², Quv=−kuv.
```

この全四 one-form coefficient は正確に
αQ=(q−dq+−q+dq−)/k である。
A-row の t0u²,t1v² は chart correction で cubic以上なので二次最低項に寄与しない。

独立 quadric pencil <F,G>→[GdF−FdG] は Gr(2,Sym² V*) 上で全域regular。
零になるなら d(F/G)=0 で独立性に反する。
従って proper image とその affine cone（零を含む）は closed GL4-invariant。
「ある order二以上zeroでこの cone条件」という incidence を P4 上で忘れる射はproper。
全一般 T36/PGL源がこの条件を満たすので whole closed T36 も満たす。
period/edge/frame の任意境界も含み、normalization-flat 仮定はない。

一方 DOUBLE の全候補 q は L1。
αQ が quadric pencil GdF−FdG から来るとすれば、J=F/G は rational first integral。
その分母 quadric は irreducible degree三 C を含めず、Cの一般点で定義される。
Cは非零 foliation の invariant hypersurface なので、dJ は C の全 tangent方向で零。
char零より J|C は定数 c、ゆえに degree二の F−cG が C を割る。
次数より F−cG=0、独立性に反する。
これは全 t0 と infinity の双方で成立する。
従って DOUBLE は T36 の全閉 incidence の外、**whole T36 に属さない**。

## 5. T38 の whole zero-surface 条件との矛盾

true completed T38 は p,t2,t3≠0 の稠密開で

```
q1=(p/2)x1x3+t2x4², q2=−(p/2)x2x4+t3x3²,
q0=p x1x2+(4t2t3/p)x3x4,
P01=−2t0q2/p, P02=2t1q1/p, P12=q0, P13=q1, P24=q2.
```

全三 q が
[x1:x2:x3:x4]=[−2t2*s³/p:2t3*t³/p:s*t²:s²*t]
で零になることは各項の相殺から従う。
x0自由の cone は二次元 integral surface で、その全域で raw P=0。
全 t0,t1 を捨てず、全 PGL像も同じ intrinsic zero dimension≥2 を満たす。
proper universal projective zero incidence の fiber dimension 上半連続性より、
この条件は whole T38 の閉包全域へ延長する。
DOUBLE の全零台は§1で正確に一次元なので **whole T38 に属さない**。

## 6. 結論の範囲

DOUBLE の有限還元が与えた actual orbit とその閉包の密点に全T36/T38非所属を得た。
その閉包の全境界点が T36/T38 と交わらないという意味ではない。
既 known24 factor/torus filters と合わせて、残る prior24 親は T37 と Ga9 cuspidal quartic。
prior23 の全監査、unknown upper bound、成分性、scheme smoothness は別問題である。
D101/M48 など既ambient arithmeticの新しい独立再認証はしていない。
