# R=3L・m二 SF00 nil：全 canonical affine／Ore 還元の別作者監査

2026-11-08、PRIMARY G。`/root/quartic_sources` が
`PLANE_TRIPLE_LINE_M2_SF00_NIL_CANONICAL_AFFINE_REDUCTION.md`
（SHA256 `1fa54ebe2f2f221fd73f19da086a697aa408cf9931014f511aaa4a7f4f65275a`）を独立に手監査した。
全 t∈C,d∈C* と full raw36/Euler6/P2九を保持した必要還元、A=0 の全 essential三認識、
A≠0 の必要 Casimir 形を受理する。**A≠0 全 G／E は未閉**である。
新 CAS、operator、matrix、minor、jet、GB、有限 parameter 標本の実行は無い。
既存 reviewer code 2bfdf012 と全 current contracts は変更していない。

## 1. fullsource の canonical lift と二つの commuting affine 場

全 first-symbol 被覆は既別作者 receipt
`PLANE_TRIPLE_LINE_M2_SQUAREFREE_CONSTANT_ROOT_SYMBOLS_SEPARATE_REVIEW.md` の nil 全二 parameter 源を入力とする。
{b,c}=3a²、mixed columns が X=N、Y=tN+dN²、N=a∂b+b∂c、N²=a∂c である。
二次係数の normalweight から higher Puv は存在しない。
raw P1 mixed columns は normal quadratic の constant base vector、P2九は base-base のみ。
従って row divergence は Ru=Rv=0、Rbase=X+Y+r1(u,v)、r1 は normal linear。
これは任意の higher coefficient を残した fullsource 恒等式である。

canonical lift Π=P+E∧R/5 は divΠ=0 で元 projective tensor を表す。
既全 canonical-lift 補題をここへ適用する条件は同じ homogeneous quadratic と五変数である。
自立的にも、E∧[Π,Π]=0 と polynomial Koszul exactness から [Π,Π]=E∧B、
BV 恒等式で div[Π,Π]=0、homogeneity から B は Euler 外積となり、strict Jacobi が零となる。
divΠ=0 と strict Jacobi から Pfaffian 一形式 i_(Π²)vol は closed、Euler weight 五の Cartan identity は
actual F の微分をその非零定数倍にする。Pfaffian kernel 恒等式から Π#dF=0。
従って canonical lift の Casimir 性を、一般の Euler 代表へ誤って課したものではない。

K=C(u,v) 上、V1=Π#du/u、V2=Π#dv/v は base-affine。
Πuv=0 と strict Jacobi により commute する。
canonical Euler correction は mixed column から (X+Y)/5 を引くので

`linV1=((4−t)N−dN²)/5, linV2=((4t−1)N+4dN²)/5`。

二係数 determinant は 3d/5、d非零から全 t で可逆である。
従って V=(4V1+V2)/3=N+p、
W=((4−t)V2−(4t−1)V1)/(3d)=N²+q。
commuting affine 場の translation 部は Nq−N²p=0、すなわち q_a=0,q_b=p_a。
A=p_a、s=(4−t)/5、r=(4t−1)/5 と置けば (V1)_a=sA,(V2)_a=rA、4s+r=3。
t=1,4,1/4 は全てこの式で覆われる。
u,v の反転は fraction-field identity 用で、これを ambient PGL 変換とみなさない。

## 2. c Γ の全域消滅と literal polynomial Ore 包含

canonical Πab の全 possible shape は
`h a²+aL+bM+cΓ+F2`、h=(1+t)/5、L,M,Γ は normal linear、F2 は normal quadratic。
strict Jac(a,b,u)/u は
`{a,(V1)b}−{b,(V1)a}−V1(Πab)=0`。

最初の項は、(V1)b=s a+translation、Πau,Πav が normal-only なので b 非依存。
第二項も (V1)a=sA が normal-only、Πbu,Πbv は a と translation のみなので b 非依存。
V1(Πab) の b 係数だけが sΓ、従って sΓ=0。
v の式も同様に rΓ=0。4s+r=3 なので全 t で Γ=0。
二つの特殊 t を generic な単位で割って消していない。

Rbase は c 非依存なので E∧R の ab entry も c 非依存、従って raw Pab の c 係数も零。
raw の他の c 依存は Pac,Pbc の c 倍 normal linear だけであり、
T=(∂cPac)∂a+(∂cPbc)∂b は global polynomial、c 非依存である。
P'=P+E∧T はこの二 c 項を同時に消し、他の full columns と P2九は全て保たれる。
よって文字通り `P'=Q(a,b,u,v)+∂c∧W0(a,b,u,v)` へ全 source が入る。

E=(E_y,c∂c) と元 projective Jacobi の c 係数を比較すれば [Q,Q]=0、
constant 係数は E_y∧[Q,W0]=0。
この whole necessary Ore 包含は rational translation の正則性に依存しない。
full Q/W の独立係数や可逆 coordinate map、G 上界が既に得られたとは呼ばない。

## 3. actual Casimir の全 A 零／非零

V(F)=W(F)=0。
A=0 では W=(a+q_c)∂c、V=(a+p_b)∂b+(b+p_c)∂c。
K[a,b,c] は domain、a+q_c と a+p_b は非零 polynomial なので F_c=F_b=0。
元 F は polynomial のため、この恒等式は u=0,v=0 の両 divisor も含めて成り立つ。
従って F は a,u,v のみの polynomial。
被約既約 quintic は essential variable 二以下では存在しないので、accepted 全 essential三・八成分定理が適用される。
この枝は全点 E 被覆と G≤23 を持つ。

A≠0 では B=p_b,C=p_c,Q=q_c、a'=a+B,b'=b+C,q'=Q−B とすると

`V=A∂a'+a'∂b'+b'∂c, W=A∂b'+(a'+q')∂c`。

I=A²c−A(a'+q')b'+a'³/3+q'a'²/2 は両場で殺される。
A が K の単位なので K[a',b',c]=K[a',b',I]、この座標で
W=A∂b'、V=A∂a'+a'∂b'。従って joint invariant ring は K[I]。
I の base-degree 三最高項 a'³/3 により、base-degree 高々五の F は I の一次式だけとなる。
original exact leading の a³ 係数は −d uv であり、
**F=−3d uv I+H5(u,v)**、H5 は homogeneous degree五 rational、という必要式を受理する。
原 note の H5 を、この座標のまま polynomial と強めない。

## 4. 独立に得られる pole-safe 追補

以下は原 note を上書きせず、本 receipt に独立の必要条件として追記する。
同じ A-linear 強化と constant-subtracted 展開は producer も別途 root へ通知している。

まず F_c=−3d uv A² は polynomial。uv は squarefree なので、各 irreducible prime ℓ の valuation で
`valℓ(uv)+2 valℓ(A)≥0`、valℓ(uv) は零又は一、ゆえ valℓ(A)≥0。
従って **A∈C[u,v]**、normal homogeneous degree 一より A=αu+βv。
これは t=4,1/4 にも有効である。
元 P0 を保つ actual normal diagonals u'=λu,v'=μv で係数は α/λ,β/μ となり、
A=u、v、u+v の全三零／非零配置へ逆を持って入る。
A=0 は前節の別枝で、d非零も保持する。

Πbu/u=s(a+B)−dA/5、Πbv/v=r(a+B)+4dA/5。
sr非零では元 mixed polynomial 性から uB と vB が polynomial、gcd(u,v)=1 より
**B は normal linear polynomial**。
t=4 では vB polynomial だけ、t=1/4 では uB polynomial だけを保持し、
この二壁の個々の rational pole を消さない。

さらに I0=I−I|a=b=c=0 と定義すると

`I0=a³/3−Aab+A²c+(B+Q)a²/2+(BQ−AC)a−AQb`。

これは直接展開による全 rational identity。
従って

`F=−3d uv I0+g5(u,v), g5=F|a=b=c=0∈C[u,v]_5`。

この **g5 は genuine binary quintic** である。
I の定数項を移さない原 H5 の polynomial 性とは区別する。
F の係数比較が強制する pole-safe 必要条件は

`uv(B+Q)∈C[u,v]_3, uv(BQ−AC)∈C[u,v]_4, uv A Q∈C[u,v]_4`。

B,C,Q 個々が全域で regular だとは主張しない。
full mixed polynomial 性は、さらに

`u(sC−dQ/5)∈C[u,v]_2, v(rC+4dQ/5)∈C[u,v]_2`

を与える。二 translation はこれら numerator と非零 3d/5 の逆から復元出来る。
これらの有限係数必要条件は次の whole affine source の設計に使えるが、
canonical divergence、strict full Q/module equations、元 F の polynomial 性の全条件に代わるものではない。

## 5. 最終 scope

全 t,d非零、全 raw36/Euler6/P2九の polynomial Ore 包含、A=0 の全 E/G23、
A≠0 の正確な necessary invariant 形を独立受理する。
A-linear の三 actual 配置、B の generic regular 性と二明示壁、g5 の polynomial 性は
全 normal divisor を保持した追加の必要構造である。
**A≠0 の存在、全 G 上界、全 E、finite full Q/W source bijection は未解決**。
この pureproof receipt は後続の選んだ点、minor 又は jet を正当化する代替ではない。
