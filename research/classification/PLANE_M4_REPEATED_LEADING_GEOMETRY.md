# rank 二 plane・m4 repeated leading quartic の全幾何入口

2026-11-07。主義務 G。実 D は既約被約 quintic、S は restriction rank 二の plane とする。
S 沿いの位数は四で、normal u,v による最初の項を

\[
 D=G_4(x,y,z;u,v)+g_5(u,v),\qquad
 G_4\in H^0(O_S(1))\otimes\operatorname{Sym}^4\mathbb C^2
\]

と書く。ここでは G4 の generic binary quartic に重根がある全枝を対象とする。
squarefree 枝の再計算はしない。
本ノートは新しい full-Jacobi/minor/GB/jet を実行せず、必要な全 leading source を幾何的に有限化する。

## 1. 任意の image-rank 一は既に八成分へ入る

leading coefficient map の image W が一次元なら G4=L(x,y,z)g4(u,v)。
従って実 D 自体が L,u,v の三変数だけに依存する。
実 D は既約 quintic なので二変数以下ではなく、exactly essential-three である。
`ESSENTIAL_THREE_QUINTIC_EIGHT_COMPONENTS.md` の全 source theorem により、
この全 Poisson 源は既受理の八つの projective23成分で被覆される。
特に全 PGL 像は次元高々23で G 非優越。

この含意は root multiplicity、R の被約性、接続の log 性や一意性を一切使わない。
generic-j/rank一だけでなく、**全 m4 image-rank 一**に適用できる。
normal quartic の任意の repeated-root 境界も落とさない。

## 2. 全 singular linear binary-quartic spaces の固定重根

有限次元線形空間 W⊂Sym4 C2 の generic member が repeated root を持つとする。
h=gcd(W) を二変数 homogeneous gcd とし W=hW' と書く。
W' の全 gcd は一である。
dim W≥2 のとき、二つの一般元 a,b∈W' は gcd(a,b)=1 と選べる。
標数零の nonconstant rational map a/b:P1→P1 の branch values は有限なので、
一般の a+t b は reduced zeros を持ち、さらに h の有限 root を避ける。
従って generic h(a+t b) の全 multiple roots は h の fixed multiple roots だけである。
これは P1 の elementary Bertini 証明であり、curve pencil の classification を仮定しない。

よって generic root partitions ごとに全 W は次の通り。

| generic partition | 全 source の必要形 |
|---|---|
| 2+2 | h に二つの fixed double roots、deg h=4、dim W=1 |
| 4 | h に fixed fourth root、dim W=1 |
| 3+1 | fixed triple root l、G4=l³(A l+B m)、dim W≤2 |
| 2+1+1 | fixed double root l、G4=l²Q、Q は generic squarefree quadratic、l∤Q |

従って 2+2 と4の **全枝**は§1の八成分/G≤23へ還元済み。
3+1 の image-rank 二と、2+1+1 の image-rank 二・三だけを残す。
3+1 では twisted GL2 Lie stabilizer が一次元あるため、接続一意性を使わない。
以下の 2+1+1 に限定する有限化と混同しない。

## 3. 2+1+1 の二つの intrinsic collision divisors

定数 normal GL2 で l=u とし、Q=A u²+2Buv+C v²、A,B,C は base 一次式とする。
Q は generically squarefree かつ u∤Q なので、

\[
 C\not\equiv0,\qquad T=AC-B^2\not\equiv0.
\]

C=0 は fixed double root と残る root の衝突、T=0 は残る二 simple roots 同士の衝突である。
既受理の associated-graded Hamiltonian 補題で G4/f は flat connection に parallel。
U=S\R 上では nonzero parallel form の root partition が変わらないため、
**両 divisor の support は R の成分だけに含まれる。**
basepoint で G4 が零となる場合も、U 上で parallel section が一点零なら恒等零となるので除かれる。
support=R は仮定せず、unused components を許す。

rank(A,B,C)=3 なら T は nondegenerate irreducible conic。
C=0 は symmetric determinant conic の tangent line である。
従って R が reduced cubic の場合、R は必ず tangent conic+line、
全 source は定数 base/normal GL により

\[
 f=z(xz-y^2),\qquad G_4=u^2(xu^2+2yuv+zv^2)
\]

へ戻る。全 nonzero 共通振幅を残す。
transverse conic+line や三直線にこの rank 三枝は存在しない。

## 4. rank 二 quadratic pencil の全二型と marked double root

rank(A,B,C)=2 なら WQ は regular symmetric2×2 matrix pencil。
det≢0 の元を一つ選び、その inverse と他の元の積を Jordan 化すれば、
全 congruence normal forms は次の二つで尽くされる。

- det の二 root が相異なる：Q=L r²+M s²、L,M independent。
- det が double root：Q=L r²+2M rs、L,M independent。

これは一般値だけの対角化ではない。第二型は symmetric Jordan pencil として保持する。
det が恒等零の pencil は Q の generic squarefree 条件に反し、ここには無い。

第一型で marked l=αr+βs とする。
αβ=0 なら交換後 l=r、全 G は r²(L r²+M s²)。
二 collision divisors は M と LM なので二つの used lines が必要である。
αβ≠0 なら diagonal rescaling で l=r+s と出来、collision C は L+M。
L,M,L+M の三相異 used lines が必要なので R は concurrent triangle である。
定数尺度と ordered-line permutation により L=x,M=−y と取れる。

第二型で l=r なら G は generic triple-root 型になり、§2の明示 3+1残余へ送る。
l が r と独立なら pencil の upper-triangular stabilizer で l=s へ戻せる。
全 G は s²(L r²+2M rs)、collision divisors は L と M²。
やはり二つの相異なる used lines が必要である。
従って conic+line 二型では rank 二枝は存在しない。

## 5. reduced reducible R の全高-rank 2+1+1 source は六型

§§3–4から、全 rank 二・三 source は次の六型に尽くされる。
前三つの二振幅表示では A,B を任意の非零定数として保持する。
三 used lines の型は collision divisor が第三 R component と一致するので比が固定され、共通非零振幅だけが残る。

| R | 全 leading generators / nonzero amplitudes |
|---|---|
| nonconcurrent xyz | A xu⁴+B yu²v² |
| nonconcurrent xyz | A xu²v²+B yuv³ |
| concurrent xy(x−y) | A xu⁴+B yu²v² |
| concurrent xy(x−y) | A xu²v²+B yuv³ |
| concurrent xy(x−y) | A (u+v)²(xu²−yv²) |
| tangent conic+line z(xz−y²) | A u²(xu²+2yuv+zv²) |

有限 line permutation、normal GL2 と本節の全非零振幅が全源を覆う。
重根の衝突や振幅零による partition 変化は、generic partition による§2の別枝として保持される。
同じ fixed connection を取る際は、全 parallel coefficient kernel を使い振幅を chosen point に固定しない。

2+1+1 quartic は三つの相異 root を持つので twisted GL2 の infinitesimal stabilizer は零。
したがって、この六 full leading normal form に parallel な full rational connection があれば一意である。
これにより nonlog/off-diagonal 接続を捨てず、将来 full30/9 source を有限個の固定 P0へ還元できる。
ただし本ノートは各 P0 の existence、全 source dimension、rank 四の非空性をまだ計算していない。
新しい finite job はこの六型を対象とする場合に別契約を作る。

## 6. 完了範囲と残余

数学的還元として全 m4 rank 一、generic22、generic4は既八成分で被覆される。
reduced reducible R の generic211高-rank全源は上の六型へ有限化された。
generic31 rank二の非一意接続、m2/m3、nonreducedRの高-rank leading、他支持は未解決。
六型の全 Poisson source、G非優越又は E帰属は未主張。
本ノートの新しい全幾何 finite reduction は独立読解を要する。

## 7. 全六型の一意接続の明示式

l=u、Q=A u²+2Buv+C v²、T=AC−B² とする。
normal covectors に作用する connection matrix を β とし、parallel convention を
d(G/f)+(βn)∂n(G/f)−tr(β)G/f=0 とする。
fixed double root を保存するため βuv=0 であり、残る三つの全係数比較から一意に

\[
\begin{aligned}
 \beta_{uu}&=\tfrac12d\log f-\tfrac14d\log T,\\
 \beta_{vv}&=\tfrac12d\log f-d\log C+\tfrac14d\log T,\\
 \beta_{vu}&=(\tfrac B2d\log T-dB)/C
\end{aligned}
\]

を得る。Euler contraction は I2、curvature convention は dβ−β∧β=0。
六型の全てで fβ の各成分は degree 二 polynomial one-form となる。
Jordan 型で T が square になっても、off-diagonal numerator の B が余分な分母を消すため、
simple-pole 接続の full source 内にある。

τ=β−(dlogf/3)I は horizontal である。
各 entry に対し iE iV vol=3fτ を満たす tracefree linear base field V を取る。
この H0Ω1(3) と sl3 の同型を完全八係数線形系で扱えば、
raw mixed brackets は {xi,na}=−Σb nb(Vab)i となる。
actual D2=0 は raw normal-normal quadratic K の全三係数を一意に定める。
全350 Jacobi/Pf0、全十五 parallel-leading coefficients と全30/9 source を検査する有限入口が得られる。
ここでまだ新算術の成功を主張しない。

fixed P0 の全源を保つ実 connected group は、三角では base diagonal 二＋normal diagonal 二＋shear六で十、
concurrent の前二型では base三＋normal diagonal二＋shear六で十一、
concurrent 第三型では base三＋normal scalar一＋shear六で十である。
tangent conic 型は base の conic+tangent stabilizer 二が symmetric-square表示により normal に持ち上がり、
normal scalar一とshear六を加えて九となる。
最後の型で base action と normal action を独立と誤算せず、実際の同時作用だけを使う。
固定 P0 の source dimension 上界 d が得られれば、各全像上界は24+d−g。
concurrent で全源の共通外積因子が得られる場合は、より弱く十分な rank四空性を直接使う。
