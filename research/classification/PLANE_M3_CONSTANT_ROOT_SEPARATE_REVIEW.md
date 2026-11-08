# plane m三・定数normal root障害の独立全源証明確認

2026-11-07、PRIMARY G、`/root/log_manuscript`。
`PLANE_M3_CONSTANT_ROOT_OBSTRUCTION.md`
（SHA256 `eb847d2c28371aa836477f713b76f5b249680389c31926f22b1794bf9994f320`）と
全root入口 `PLANE_REDUCED_CUBIC_M3_ROOT_GEOMETRY.md`
（`6c665327f312bdc371b9c099ac39b2b68fbb70e3e1467bd127a09ced7554ea24`）を別作者として全読解した。
以下のwhole必要条件・明示空性・first-connection必要還元を **独立受理** する。
新しいmatrix、CAS、anchor、jet又はGBは一切実行していない。
全m三又は全G/Eの完成、残connectionの実Poisson liftの存在は結論しない。

## 1. 全36/Euler6と全9から直接得る係数

baseをa,b,c、normalをu,v、vol=da∧db∧dc、base EulerをEとする。
Pi_S=P_f は `{a,b}=f_c,{b,c}=f_a,{c,a}=f_b`。
constant normal rootをuに取り、全tracefree mixed gaugeを戻せば

```
P0_base=P_f,
P0^{i,u}=u X_i, P0^{i,v}=u Z_i+v Y_i,
P0^{u,v}=k20 u²+k11 uv
```

である。root保存はu-columnのv項を零にする。
normal bracketのv²係数はD2のv²係数へ3f倍で寄与し、G2=0から零。
Z、k20、k11を消したchosen sliceではない。

全P1は、base-baseの三pair×base一次三×normal二の18係数と、
mixed六pair×normal二次三の18係数、合計36係数からなる。
Euler gauge六は `E_total∧(u A+v B)`、A,Bは任意constant base vectorsである。
そのi,u成分は `−u(u A_i+v B_i)` なので、bad coefficient
`P1^{i,u}|u0=v²C_i` **はこの六方向で変わらない**。
従って忠実P1三十に移ってもCを失わず、raw代表の選択に依存した架空の障害ではない。
全P2はbase-base三pair×normal二次三の九係数で、normal columnを持たない。

u零上で、v³へ寄与し得る最低項は

```
P_base=P_f+O(v),  P^{i,u}=v²C_i,
P^{i,v}=vY_i+O(v²),  P^{u,v}=0.
```

Euler-signed actual Dを、座標順(a,b,c,u,v)の `E_total∧P∧P/2` とする。
Eulerのv項によるPfは
`v*(f_c v²C_c+f_b v²C_b+f_a v²C_a)=v³ C(f)`。
base Eulerと二mixed columnの三項は `−v³ det(E,C,Y)`。
P1 base-base、他のmixed P1、全P2はこの係数でnormal degreeが四以上、
P0のZ/Kはuを含むので寄与しない。
よって全sourceで厳密に

```
[D3]v³=C(f)−det(E,C,Y).                         (A)
```

C零なら全u-rowはu倍となり、u=0上のprojective tensorは三次元tangent空間内でrank二以下。
従ってactualDはuで割れ、被約既約quinticに反する。全対象でC非零が必要である。

## 2. connection characterとfactor三の独立確認

c一chartを取り、intrinsic normal coordinateはv/cである。
quotient bracketを直接微分すると、v係数の二Hamiltonian actionsは

```
{a,v}=v(Ya−aYc+f_b)+u(...),
{b,v}=v(Yb−bYc−f_a)+u(...),
{a,b}=3f.
```

従ってcomplement characterは

```
q=[(f_a−Yb+bYc) da+(f_b+Ya−aYc) db]/(3f).
```

homogeneous表示は
`3f q=df−i_E i_Y vol`、q(E)=1である。
実際 `(i_E i_Y vol)(C)=det(E,C,Y)` なので(A)は

```
[D3]v³=3f q(C).
```

normal root uによりG3|u0=0、f非零なので、**全constant-root源で q(C)=0、C≠0**。
ここまではRの被約性、root multiplicity、他rootのsplit性を仮定していない。
projective quotientのEuler補正を落としたraw mixed fieldの推測でも、
追加strict Jacobiを課す議論でもない。

## 3. reduced Rに限ったactive-prime比較

R被約の場合に限り、triangular flatnessはdiagonal qをclosedにし、元fβのpolynomial性からpoleは高々単極。
各primeでのclosednessからresidueはconstant、全residueを引いた差はP2上holomorphic一形式となって零。
従って `q=Σρ_i dlog f_i`、`Σdeg(f_i)ρ_i=1` で尽くされる。
q(C)=0をprime f_iで比較すれば、active ρ_i非零では `f_i|C(f_i)`。
右辺は次数一小さいからC(f_i)=0である。

不可約conic Qはrank三なので、C(Q)=0はC零を強制する。
従ってR=LQでは **q=dlogL、C(L)=0** だけが許される。
independent triangleでは三active residuesはC零を強制し、少なくとも一residueが零でなければならない。
concurrentでは全lineを殺す共通translation方向があるため、この部分を空性へ延長しない。
**非被約Rにはこのsimple-pole/log-residue分類を適用しない。** 使えるのは前節のq(C)=0とC非零までである。

## 4. 全root degreesと21の全零residue補集合

全source入口のUFD/二gradingを確認した。
generic21はmultiplicityの違いで両rootがrational、primitive度(d,e)のbudget `2d+e≤2` から
00/01/02/10の全四flagsである。
tripleは3d≤2なので必ずconstant root。
SF splitは000/001/002/011の全四flagsであり、SF nonsplitをこの表で代用しない。
nonconcurrentのπ1はabelianなのでSF imageは1/C2/C3、concurrentのS3は明示残余。
splitting cover上のbudget `ΣL_i·p*H≤2 deg p` はprimitive contentとnef性のみであり、
cover compactificationのflatness又はPicard零を仮定していない。

21のrow frame Mから `dM+Mβ=diag(θd,θs)M` を代入すると
`θd=dlog(f/(h detM))`、`θs(E)=e+1`、θsは全closed character。
moving degree一のoff-diagonal `(bη−db)/L` では、primitive bと単極性からηのL-residueが零。
moving degree二でdetM=L²の場合はtangentだけでなく全normal derivativeを比較する。
conicでA=cQ+dL²、triangleでA=cM^nN^(2−n)+dL²、n=0,1,2が全てとなる。
これらの既全式を戻してcomplement q又はpのactive residuesを照合した結果は次の通り。

| R / partition | 全必要結論 |
|---|---|
| conic+line / 21 / 01,02,10 | conic-active characterにより全空 |
| conic+line / 21 / 00 | h=Q、p=q=dlogLだけ。h=L²は空 |
| triangle / 21 / 02,10 | 三active residuesにより全空 |
| triangle / 21 / 00 | h=LiLj、p=dlogLk、qの全三residue-zero affine lines |
| triangle / 21 / 01、h=L | η=ρdlogM+(1−ρ)dlogNのρ=−1又は2 |
| triangle / 21 / 01、h=M | 同ρ=0又は2。h=Nは交換した両値 |

例えばconic01は `q=(3/2)dlogQ−2dlogL`、02は
`3dlogL−dlogQ` 又は `2dlogQ−3dlogL`、10のconstant simple rootのcomplementは
`dlogQ−dlogL`。いずれもconic-activeである。
triangle02の平方determinantでは
`q=−3dlogLi+(1+n)dlogLj+(3−n)dlogLk`、n=0,1,2で全三active。
01の許すfinite periodsでも全linear bとresidue零で増えるC方向を保持し、full sourceを一点へ置換しない。

## 5. split SF、C2 fixed-constant、tripleの必要入口

split SFの000ではβ=dlog(f/h)I。
conicはh=Qだけ、triangleはh=LiLjだけ残る。
001はconicでは存在せず、triangleではhが二moving coefficientsのunused第三lineの場合だけ、
p=dlogb、q=dlogaとなる。h=a/bは一方に三active residuesを与える。
002は、conicではQ-active、triangleでは少なくとも一平方quadraticにより三activeとなって全空。
011のb1=b2=Lも同じcomplement dlog(f/L²)で全空、distinct b1,b2のtriangleは
determinantが残第三lineの平方かつ二lineのideal内というbasepoint矛盾で空。
concurrent及びSF nonsplitの全complementsは保持する。

C2でfixed rootがconstantなら `G3=u(Au²+2Buv+Cv²)`、q=dlog(f/C)。
conicではC=const Q。discriminant T=AC−B²の全degree四支持はL⁴/L²Q/Q²であり、
nonsplit性は平方二型を除く。従ってT=τL²Q、τ非零。
Q|B²からB=bQ、constant shearでB零となり、全leadingは
`G3=u(τL²u²+Qv²)`、`β=diag((1/2)dlogQ,dlogL)`。
SF rigidityによりoff-diagonalも零である。amplitudeをhigher sourceから消す推論は用いない。
moving fixed-root C2とconnected C3/S3は別残余。

tripleでは `q=2p−dlog(f/h)`、したがってconicの全constant rootは
h=Qでp=q=dlogL、h=L²でp=(1/2)dlogQ、q=dlogL。
全nonlog off-diagonalをr=α_Z、Z tracefree linearとして保持する。
row-connection conventionのflatnessは `dr+(p−q)∧r=0`。
直接Cartan計算は

```
dα_Z=Z(f) i_Evol/f²,
η∧α_Z=−η(Z) i_Evol/f  (η(E)=0)
```

を与え、従って全式は **(p+dlogh)(Z)=0**。
h=QではZ(f)=0でrはdlogQ−2dlogLの任意constant倍。
h=L²では `4Q Z(L)+L Z(Q)=0`、互いに素性から
Z(L)=aL、Z(Q)=−4aQ。Q非退化とtrace Z零からa零。
ゆえにrは(1/2)dlogQ−dlogLの任意constant倍で、constant normal unipotentにより消せる。
前者のnilpotent extensionの零/非零を両方残す。
これは全八係数rからの必要還元であり、初めからlog rだけを選んでいない。

## 6. 正確な受理範囲

全source係数障害と上表の明示空性・first connection必要還元を受理する。
残scalar/diagonal/nilpotent first symbolsがfull P0を持つか、全higher liftのG上界又はE帰属は別義務。
concurrent、nonsplit SFの残covers、三線triple等を閉じていない。
非被約Rの既完了packetへ今回の単極residue比較を無断で適用しない。
本receiptは独立手証明であり、計画gate又は機械演算が幾何的主張を自動証明したという意味ではない。
