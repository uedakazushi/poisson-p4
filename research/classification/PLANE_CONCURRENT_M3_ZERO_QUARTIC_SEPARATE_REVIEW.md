# Concurrent m三・K4零全Ore源の独立被覆確認

2026-11-08、PRIMARY G/E、別作者 `/root/quartic_review`。
constant primitive direction、primitive invertible linear direction の SS/Jordan、
constant-proportional pair と zero pair を含む **whole K4零 source** を独立読解した。
指定 rank 四・actualD 被約既約 quintic の全点は essential-three 既受理八成分へ入る。
その global 像次元は高々23。新成分、exact orbitdimension、ambient componenthood は主張しない。
新 CAS、数値実行、operator/module、source probe、minor、jet 又は GB は無い。

## 1. 全源、共通方向の次数と補集合

対象は rank 二 invariant plane、restriction cubic が三被約 concurrent lines、normal 位数三、
全 first-cover の Ore 部分の全 higher。
既独立確認済み strict/Euler bridge により

    P=Q(y)+∂c∧W(y), y=(a,b,u,v),
    Qab=aA1+bB1+C2, Qan=U2,Qbn=V2,Quv=0,
    [U2,V2]=A1U2+B1V2, E_y∧L_WQ=0。

全 raw36/Euler6 P1、全九 P2 の affine particular、全 C2/W を保持する。
K4=PfQ=0 なので actualD は c 非依存。
この receipt は元 non-Ore first profiles を含むと仮定しない。

非零 pair は UFD で U=P_(2−e)R_e,V=Q_(2−e)R_e、primitive R_e、e=0,1,2 と書く。
P/Q が constant-proportional でなければ e≤1。
e0 の binary quadratic pencil は common-root/basepoint-free の二類と、
marked constant direction の二位置ずつで C_tan/C_trans/B_axis/B_other に尽くされる。
direction を無視して normal frame を潰していない。
e1 は rank 一 matrix なら共通 linear factor があり primitive 性に反する。
従って invertible linear direction の SS(t!=0) 又は Jordan に尽くされる。
SS の zero eigenvalue は C_trans、nilpotent degeneration は C_tan へ戻す。
全 normal/base frames は actual constant GL blocks、全 inverse は通常の matrix inverse。

central の除去は actual positive base shift の四列を使い、rank-drop source を残す。
この小さい正常化を full W の equation と取り違えない。
変換後も全 W が quadratic であること、元 normalweight/actualD 条件を戻すことを保持する。

## 2. constant-proportional pair と zero pair

U=λR,V=μR、(λ,μ)!=0 なら strict Jacobi から λA1+μB1=0。
v0=λ∂a+μ∂b と置くと v0∧Q=0,L_v0Q=0。
actualD は符号を除き E_y∧W∧Q の係数であり c 非依存。
v0 微分では v0(E_y)=v0 の項が消える。
module を v0 で wedge した式から残りの E_y∧[v0,W]∧Q も零。
従って D は c と v0 の二独立 constant directions に依存しない。
非零既約 quintic は essential 二以下ではないので essential 三。
既受理全 essential-three theorem により全点が八成分へ入り、global 上界23をそのまま使う。
plane/frame/v0 選択の次元を重複加算しない。
この constant-factor lemma は既独立 receipt の全源証明にも一致する。

U=V0 では Q=F2∂a∧∂b、actualD は F2 times cubic、F2零なら D零。
全 C2/A1/B1 を含め rank 四既約 quintic の対象外。

## 3. constant direction：unimodular 部分と全 C_tan period

B_axis/B_other 全 C2、C_tan/C_trans の t0 全 C2 は
normal linear Casimir L と divQ0 を持つ。
chart L1 では全 module が L_Ypi=0、div pi0 から pi=Jac(F)、
F=aQ2(r)−bP2(r)+C3(r)、degree d≤3。
全 central は C3' として保持される。
module は d(YF)=(divY)dF を与え、C(F,b,r) 又は C(F,a,r) の正確な field 表示から YF∈C(F)。
polynomial YF の rational denominator は Bezout substitution で消えるので YF=P(F)、P polynomial。
degree五の homogenization は L^5P(F_d/L^d)。d3 なら L²、d2 なら L 因子、
d1 は二 linear forms の binary quintic。zero complements も既約 rank 四対象外。
この部分は先行 nonzero-K receipt §7 でも独立確認済みである。

C_tan t!=0 では shift map rank三なので全 C2 を actual に零へ戻せる。
chart u1、r=v、x=b−ar は chart polynomial automorphism、逆 b=x+ar。
元 pi=∂a∧V、V=∂r−t x∂x を直接確認した。
全 polynomial Poisson field Y の H=Yx+t xYr は H_a0,VH=−tH。
H=Σ h_j(r)x^j なら h'_j=t(j−1)h_j。
t!=0 の全 polynomials で j1/constant のみなので H=λx。
actualD=λu³(bu−av) 又は零。全 period・全 W で proper factor。
chart polynomial automorphism を ambient PGL と呼ばない。

## 4. C_trans：全複素 period と三 central walls

chart v1、r=u、x=a/r−b、k=t+2 の dense open では
pi=∂b∧V、V=r∂r−(kx+C2(r)/r)∂x。
H=rYx+(kx+C2/r)Yr に対し H_b0,VH=(1−k)H。
元 chart volume と新 volume の比は r なので actualD0=±rH、
V D0=(2−k)D0。
z=a−br=rx へ戻すと polynomial ring C[b,z,r] 上で

    V=r∂r+(1−k)z∂z−C2(r)∂z, D0_b=0。

normal 位数≥3 は base-degree≤2、従って z-degree≤2。
C2零なら可能 monomials は r^(2−k), r z, r^k z² の全て。
指数が非負整数の場合だけ許すので全複素 k を扱う。
u=r 因子を避けるには k0又は2のみ。
k0 は元 a,b,r の ordinary degree≤4 なので quintic homogenization に v 因子。
k2 の r²z² は degree六の b²r⁴ を他項で消せず、global degree五からその係数零、残りも v 因子。
他の全 k では u 因子。

nonzero central residual は (k,C2)=(-1,r²),(0,r),(1,1) の三 whole walls。
S=r∂r+(1−k)z∂z と N=−C2(r)∂z は可換、N は各 polynomial 上 locally nilpotent。
S の weight space ごとに nonzero scalar+N が可逆なので、必要 eigenfunction D0 は所定 weight かつ ND0=0。
domain 性から D0_z0、D0_b0 より D0 は rだけ、全三場合は r³/r²/r の倍数。
全 central zero/nonzero と全 W が proper factor 又は zero で排除された。

## 5. SS：全 module、全 period と global D の完全 monomial budget

全 SS(t!=0,s任意) と全 C2 を使う。
chart u1、x=a,h=b−ra,r=v の literal dehomogenization は

    pi=∂x∧V, V=(μh+C2(1,r))∂h+νr∂r,
    μ=−s−2, ν=t−1, Z=∂x。

元 full module は L_Ypi+Y∧∂x=0。
J=det(V,Ynormal) と置くと J_x0、VJ=(divV−1)J=(μ+ν−1)J。
vertical twist の −1 を落として普通の Poisson vector field と扱っていない。
actualD|u1=±J、全元 tensor の Jacobi から得る necessary identity である。

三 central walls は μ=jν,j=0,1,2、residue C2=γr^j。
t!=1,γ!=0 では diagonal part と γr^j∂h が可換 locally nilpotent。
weight 分解から J_h0、J_x0、従って homogeneous D は binary(u,v) quintic で可約。
γ0 は次の全 diagonal source へ戻す。
t1,s!=-2 は ν0、μ!=0、J の h-degree j≤2。
μ(j−1)=−1 は j1 を許さず、j0 は binary、j2 は (bu−av)² times binary linear。
common center(t1,s−2)は全 C2 三係数を残して Q=−E_y∧D0+C2∂a∧∂b、actualD は C2 times cubic 又は零。

C2零,t!=1 では global quintic の全項は w=bu−av として

    D=Σ c_jk w^j v^k u^(5−2j−k),
    j=0,1,2, k≥0, 2j+k≤5。

negative u-exponent を勝手に捨てていない。
最悪 pole の u0 coefficient は異なる a^j v^(j+k) の和なので cancellation 不可。
このため全項が指定 budget 内へ入る。
normal 位数正確三は j2 の非零項を要求し、k0なら s=−t、k1なら s=−1。
t!=1 では両方同時にはない。

s=−t、q=(t−2)/(t−1) では全 possible terms は
u w²、u^(3−q)v^q w、u^(5−2q)v^(2q)。各指数の非負整数条件を保持する。
proper w factor を避けるには第三項が必要。
proper u factor を避ける唯一候補は 2q=5。q3 は第三項が許されず w factor。
従って (t,s)=(1/3,−1/3)、D=A u w²+B v⁵、AB!=0 のみ。
s=−1 は actual u/v,a/b swap と tensor 1/t により (t,s)→(1/t,s/t) なので、
唯一の対称候補は(3,−1)。
これは全 period の帰結であり、f の roots を a,b に仮定していない。
f の a³/b³、zero residue、他 marked cubic の補集合を除く入口も不要である。

最後の候補では μ=−5/3,ν=−2/3、J=A h²+B r⁵。
x=h0 の determinant は J=−νrYh なので Yh=(3/2)B r⁴。
しかし元の global quadratic W の exact chart formula から同じ locus で

    Yh=Wb(0,0,1,r)−r Wa(0,0,1,r)

は r-degree≤3。
B!=0 と矛盾、B0 は proper factor。swap候補も同じ。
全 SS source が空であり、未確認の f-root/coefficient classification を使わない。

## 6. Jordan：全 period と critical central residue

全 Jordan Q と全 C2 は D0∧S と producer の full linear S で正確に表される。
chart u1、x=a,h=b−ra,r=v では

    pi=∂x∧H, H=(μh+C2(1,r))∂h+∂r, μ=−s−2。

radial decomposition の Z=−∂x と W の homogeneity を戻すと、
horizontal module は [Y,pi]+Y∧∂x=0。
SS と同じ determinant F は F_x0、H(F)=(divH−1)F。
私は別 chart v1 からも同式を導出し、符号・critical parameters が一致することを確認した。

s!=-2 で C2 は rank三 shift により零。
F=Σ h^j f_j(r) の各 coefficient は f'_j=(μ−1−μj)f_j。
polynomial は f_j constant、μ(j−1)=−1 の一つの j だけ。
global degree五は D=λu^(5−2j)(bu−av)^j、regularity から j≤2。
j1 は equation に不可能、j0/j2 は proper factor、zero F も対象外。
μ0 なら F_r=−F の polynomial solution は零。

s−2 の full residual C2=γv²、γ0/nonzero を残すと H=∂r+γr²∂h。
polynomial inverse を持つ chart change h'=h−γr³/3 により H=∂r。
F_r=−F は零なのでこの全 critical magnetic source も actualD零。
chart ODE 用変換を ambient PGL normalizer と扱っていない。
全 Jordan periods/central/W が empty。

## 7. 受理した全 K0 結論と他枝

nonproportional U/V の全 primitive constant direction と invertible linear direction は empty。
U=V0 も empty。唯一残る constant-proportional source は essential 三の既受理八成分へ全点被覆。
従って指定 **whole K4零 Ore source は pointwise E8、global 像次元≤23**。
この E8 判定を K4非零の二軌道へ拡張しない。
先行独立 receipt の nonzero-K 全像≤24と合わせれば whole strict Ore 源の G上界24となる。
non-Ore first profiles は別義務、他 normal 位数・surfaces・全分類完成も本 receipt の主張ではない。
