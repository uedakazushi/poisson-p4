# R=3L・SF00 nil 全(t,d)：canonical affine対から全Ore包含へ

2026-11-08、PRIMARY G、`/root/quartic_review` の新しい純手証明。
**独立review前。** 新 CAS、matrix、operator、minor、jet 又は GB を実行していない。
全 first-symbol 正常化は既独立確認済みであるが、本稿の新必要還元と既約性の分岐は別の確認を要する。
G/E 全閉鎖は結論しない。

## 1. 対象と full normal source

座標(a,b,c,u,v)、rank二plane u=v=0、f=a³、normal位数二SF00の nil whole first symbols を

    {b,c}=3a²,
    {b,u}=a*u, {c,u}=b*u,
    {b,v}=t*a*v, {c,v}=(t*b+d*a)*v,
    other brackets=0,  t∈C, d∈C*.

とする。actual leading は −d*a³*u*v。
全 raw36 P1/Euler6、全9 P2、全 higher parameters を保持する。
quadratic polynomial の normalweight 支持から full Puv=0。
raw P1 の mixed normal columns は normal quadratic の constant 倍、P2 は base-base のみなので、
全 row divergence は

    Ru=Rv=0,
    Rbase=(X+Y)+r1(u,v),
    X=a∂b+b∂c=N,  Y=tN+dN²,  N²=a∂c.

r1 は normal linear の三vector。
P0 の tracefree mixed fields と全 higher の支持を使った恒等式で、generic parameter の推測ではない。

canonical strict lift Π=P+E∧R/5 は同じ projective tensor で divergence zero。
projective Jacobi と BV identity により strict Poisson であり、actual quintic F の Casimir lift である。
Πuv=0 は fullsource で保たれる。

## 2. 全(t,d)で二つの affine Hamiltonian fields を分離する

K=C(u,v) 上で V1=Π#du/u、V2=Π#dv/v を base(a,b,c) の vector fields と見る。
両者は u,v を殺し、strict Jacobi と Πuv=0 より commute。
full mixed columns と R の形から base に affine linear であり、linear parts は正確に

    lin V1=((4−t)N−dN²)/5,
    lin V2=((4t−1)N+4dN²)/5.

determinant は 3d/5!=0。従って **全 t、全 d!=0** で

    V=(4V1+V2)/3=N+p,
    W=((4−t)V2−(4t−1)V1)/(3d)=N²+q

を作れる。p,q は K 上の constant translations、normal homogeneous degree 一。
commuting condition は Nq−N²p=0、すなわち

    q_a=0, q_b=p_a.

A=p_a、s=(4−t)/5、r=(4t−1)/5 とすると、元 fields の a-components は

    (V1)_a=sA, (V2)_a=rA.

s と r は同時に零にならない。t=1,4,1/4 も同じ公式内で保持される。
u/v の反転は fraction-field identity を読むためだけで、source parameters の局所化ではない。

## 3. c依存の唯一のEuler不変項が全域で零

canonical ab coefficient の full possible shape は

    Πab=h*a²+aL(u,v)+bM(u,v)+cΓ(u,v)+F2(u,v),
    h=(1+t)/5,

L,M,Γ は normal linear、F2 は normal quadratic。
V1_b=s*a+(translation)、V1_c=s*b−(d/5)*a+(translation)。
strict Jac(a,b,u)/u は正確に

    {a,(V1)_b}−{b,(V1)_a}−V1(Πab)=0.

最初の二項は b に依存しない。
V1(Πab) の b coefficient は sΓ だけなので sΓ=0。
同様に v equation は rΓ=0。
ゆえに全 t で Γ=0。

Rbase の上記形より canonical Euler correction の ab entry は c非依存。
従って raw Pab 自体の c coefficient も零である。
raw P の残 c依存は Pac/Pbc における c 倍 normal linear terms のみ。
T=(Pac,c)∂a+(Pbc,c)∂b とすると P'=P+E∧T はそれらを同時に消す。
P2 全9は最初から c非依存で、affine particular もそのまま保持する。

従って全 source は文字通り

    P'=Q(a,b,u,v)+∂c∧W0(a,b,u,v)

へ入る。これは global polynomial Euler gauge であり、有理 field 変換を ambient PGL とみなす操作ではない。
projective Jacobi の c coefficient と constant coefficient は

    [Q,Q]=0, E_(a,b,u,v)∧[Q,W0]=0

を与える。ここまでに有限 parameter sample 又は chosen higher slice は無い。

## 4. actual Casimirの全分岐：A=0は既全essential三へ

V(F)=W(F)=0。
A=0 なら W=(a+q_c)∂c、V=(a+p_b)∂b+(b+p_c)∂c。
K[a,b,c] は domain であり、a+q_c と a+p_b は非零なので F_c=F_b=0。
元 F は polynomial だから全 coordinate divisor 上もこの恒等式を満たす。
従って actualF は(a,u,v)だけに依存し、対象既約 quintic は accepted essential-three 八成分へ全点被覆される。

A!=0 では a'=a+p_b、b'=b+p_c、q'=q_c−p_b とする。逆は同じ constants を引く affine K-map。

    V=A∂a'+a'∂b'+b'∂c,
    W=A∂b'+(a'+q')∂c.

joint invariant を

    I=A²c−A(a'+q')b'+a'³/3+q'a'²/2

と置く。直接 V(I)=W(I)=0。
A unit により K[a',b',c]=K[a',b',I]、W=A∂b'、V=A∂a'+a'∂b' なので joint invariant ring は K[I]。
base degree 高々五と I の degree三から、F は I の一次式だけである。
original normal leading −d*a³*u*v と degree一 translations により

    F=−3d*u*v*I+H5(u,v).

H5 は K 内の homogeneous degree五 rational function。
I の coefficients 自体が u/v に pole を持ち得るため、H5 の polynomial 性をここで仮定しない。
全 F の polynomial 性、full strict Q/module equations はなお必要である。
この式から A!=0 源の全存在、全 E 帰属、G 上界を結論したとは呼ばない。

## 5. 次の最弱義務

新しく証明したのは全(t,d)の literal c-independent/Ore 包含と、A=0の既全essential三認識、
A!=0の正確な cubic-invariant 必要形である。
全 source parameters・全9 P2・tの全特殊値を保持する。
今後算術を行う場合は、先に full Q/W source、A零/nonzero、全 normal pole/parameter walls と成功時の含意を固定する。
追加 minor/jet の系列、d=1 又は t の数値 sample はこの whole lemma の代替にならない。
