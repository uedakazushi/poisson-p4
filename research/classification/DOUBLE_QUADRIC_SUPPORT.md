# 二重二次因子の全支持：滑らかな場合と階数4の二配置

2026-11-07。分類義務 A/E。実際の退化因子が
\(D=H+2\mathcal Q\) である一般階数4の Poisson 構造を扱う。
以下は固定した代表の接空間計算に留まらず、指定した因子配置を持つ全点の命題である。

* \(\operatorname{rank}\mathcal Q=5\) なら H は接超平面であり、全点は既知の
  `DOUBLE_SMOOTH_QUADRIC` の22次元軌道に属する。
* \(\operatorname{rank}\mathcal Q=4\) で H が頂点を避ける場合、該当する Poisson 構造は存在しない。
* \(\operatorname{rank}\mathcal Q=4\) で H が頂点を含み、商の滑らかな二次曲面に非接触の場合、
  全点は定数非退化・一次部分零の多項式 Poisson 代表を持つ。
  従って既に証明された β=0 層の被覆へ帰着する。
* \(\operatorname{rank}\mathcal Q=4\) の接触配置も、全点が既知の MPS40 成分の和に入る。

非接触配置の被覆は `BETA_ZERO_COVERAGE.md` の22個の Q45 成分定理と、その全閉包被覆を依存に含む。
rank4 の接触配置は第7節で、上記の実一変数族の閉包と一つの追加軌道へ全点を還元する。
一変数族は独立した T18 認識に接続し、追加軌道は実 Poisson pencil と SNC 軌道極限で被覆する。
rank3 以下は本ノートの被覆結論に含めない。

## 1. 逆二形式への全点還元

\(F=Lq^2\) を退化五次式とし、\(\Pi\) を divergence-free 二次斉次 lift とする。
標準恒等式 \(\Pi^\sharp dF=0\) と \(\gcd(L,q)=1\) より
\(\Pi^\sharp dL=LX\)、\(X(L)=0\) となる線形ベクトル場 X が存在する。
\(L=h\) として \(P=\Pi-E\wedge X\) と置けば、h は P の Casimir であり、
P は同じ射影 bivector を表す。\(E\wedge[P,P]=0\) の \(\partial_h\) 係数から
\([P,P]=0\) が従う。従って h=1 上の P は次数高々二の真の affine Poisson bivector である。

\(Q=q|_{h=1}\) とする。\(\operatorname{Pf}(P)=\kappa Q^2\)、\(\kappa\ne0\) なので、
逆 symplectic 形式は、係数次数高々二の二形式 α によって
\(\omega=\alpha/Q^2\) と書ける。非退化性と閉性は

\[
 Qd\alpha-2dQ\wedge\alpha=0,\qquad
 \alpha^2/2=kQ^2\,\mathrm{vol},\quad k\ne0. \tag{1}
\]

この還元では affine chart 上の三次尾部を勝手に捨てていない。
不変超平面を用いた実際の鉛直 gauge によって二次代表を得ている。
Q=0 上で逆形式を用いなくても、得た Jacobi 恒等式と族は多項式として全域へ戻る。

## 2. 五つの配置を含む全90係数の線形補題

滑らかな projective quadric の二配置、および rank4 quadric の三配置の affine 標準形は次である。

| projective 二次式 | H の位置 | Q |
|---|---|---|
| rank5 | 非接触 | \(1+u^2+v^2+z^2+w^2\) |
| rank5 | 接触 | \(w+u^2+v^2+z^2\) |
| rank4 | 頂点を避ける | \(u^2+v^2+z^2+w^2\) |
| rank4 | 頂点を含み非接触 | \(1+u^2+v^2+z^2\) |
| rank4 | 頂点を含み接触 | \(w+u^2+v^2\) |

五つ全てで (1) の線形閉性条件の解はちょうど

\[
 \alpha=Qd\beta-dQ\wedge\beta, \qquad \deg\beta\le1. \tag{2}
\]

全90係数の線形写像 D は140×90、affine-linear 一形式20係数からの写像 B は90×20である。
新しい外微分による構成で \(DB=0\) を QQ 上で厳密に検査し、D の階数71、B の階数19を
QQ 行簡約と非零有理小行列式で再確認した。従って \(\ker D=\operatorname{im}B\)。
計算はこの五つの標準形の全空間を検査するもので、Poisson パラメータの標本ではない。

実装 `code/classification_double_quadric.py` の `linear`、証明書
`reports/CLASSIFICATION_DOUBLE_QUADRIC_FORMS_20261107.json` が行列の hash・minor・階数を保持する。
新実行は約2.35秒、55 CPU秒・1536MiB以内。旧整数係数構成と新 SymPy 外微分構成は独立である。
初回のlinear/rank3/family報告が参照するproducerの全bytesは
`reports/CLASSIFICATION_DOUBLE_QUADRIC_SOURCE_20261107.json`に保存した。
そのsource SHA256は`2b511b2fdb744cc004c76faa5dc9c3fa4b79c2f51c9c73601f43dc651d25b2cd`である。
現在の共通runnerには後続のanisotropic modeが追加されている。旧報告と旧producerは変更しない。

以下 \(\Omega=d\beta\)、\(\Omega^2/2=A\,\mathrm{vol}\) と書く。(1),(2) は

\[
 \Omega\wedge dQ\wedge\beta=(A-k)Q\,\mathrm{vol} \tag{3}
\]

を与える。\(\Omega=0\) は k≠0 に反する。Ω が symplectic のとき
\(\iota_Z\Omega=\beta\) により affine-linear ベクトル場 \(Z=Ly+v_0\) を定めると

\[
 \mathcal L_Z\Omega=\Omega,\quad \operatorname{tr}L=2,\quad
 Z(Q)=cQ,\qquad c=1-k/A\ne1. \tag{4}
\]

## 3. rank5 の全点認識

この場合は既存 `legacy/v03/proofs/smooth_double_quadric_theorem.txt` 第3–7節の証明を
全文点検し、第2節の有限補題を上記の独立計算で再確認した。
Ω が rank2 なら (3) の左辺は二つの一次式の積の差なので、斉次化した Q の rank は高々4となり矛盾。
従って Ω は symplectic。Z の h 固定斉次化と非退化な五変数 q から
\(2\operatorname{tr}L=5c\)、従って c=4/5。
非接触標準形の臨界点では Q=1、ZQ=0 となり矛盾する。

接触標準形では Hessian H の rank は3、\(S=\Omega^{-1}H\) は
\([L,S]=S/5\) を満たす。全 \(\operatorname{tr}S^m=0\) より S は nilpotent、rank3なので
一つの長さ4の Jordan block を持つ。\(\ker S\) 上の L の固有値4/5から、その filtration の
固有値は順に4/5,3/5,2/5,1/5となる。L は可逆で対角化可能である。
Z の唯一の零点へ平行移動し、固有座標で

\[
 \Omega=da\wedge dd+db\wedge dc,\quad Q=d+b^2+ac,\quad
 Z=(a\partial_a+2b\partial_b+3c\partial_c+4d\partial_d)/5
\]

へ規格化できる。従って全点が既存の一つの軌道に属する。
この軌道の全空間成分性・次元22は `ADDITIONAL_RANK4_EARLY_COMPONENT_AUDIT.md` の
全150係数の接空間と軌道次元の挟み撃ちによる。

## 4. rank4・H が頂点を避ける配置は空

ここでは Q は四変数の非退化な斉次二次式。
Ω が symplectic なら (4) の二次係数比較から
\(L^TH+HL=cH\)。H が可逆なので \(2\operatorname{tr}L=4c\) であり c=1となって矛盾する。

Ω が rank2 の場合、線形座標を変えて \(\Omega=du\wedge dv\) とする。
\(X=\beta_w\partial_z-\beta_z\partial_w\) と置くと、(3) は

\[
 X(Q)=-kQ.
\]

一方 \(\operatorname{div}X=\partial_z\beta_w-\partial_w\beta_z
=\Omega_{zw}=0\)。非退化 Hessian に対する二次係数比較から
\(2\operatorname{tr}(DX)=-4k\)、従って k=0となり矛盾する。
X の定数項の有無にこの trace 議論は依らない。Ω の全ての可能な階数を排除したので全配置が空である。

## 5. rank4・H が頂点を含む非接触配置の全被覆

Q=1+u²+v²+z²、臨界集合は affine 直線 C=(u=v=z=0)、そこで Q=1。
Ω が rank2 なら前節の XQ=−kQ を C 上で評価して k=0となる。
従って Ω は symplectic であり、(4) を C 上で評価して c=0となる。

\(ZQ=0\) より Z は C を保ち、L は Hessian の一方向 radical R を保つ。
三次元商上の L は Q の直交 Lie algebra に属するので trace0。
\(\operatorname{tr}L=2\) より R 上の固有値は2であり、C 上の affine ベクトル場は唯一の零点を持つ。
この点へ平行移動すると Z の定数項は零となる。
Q は引き続き定数+二次であり、\(\beta=\iota_Z\Omega\) は一次斉次である。
従って (2) の α、さらに逆 bivector P は定数+二次のみを持つ。
臨界点の Q=1、α=Ω から定数項は非退化である。

線形座標変換と全体の尺度で定数項を標準 J に合わせれば
\(P=J+Q'\)、かつ次数別 Jacobi から \([J,Q']=0\)、\([Q',Q']=0\)。
これは β=0 層そのものであり、`BETA_ZERO_COVERAGE.md` により既知の閉成分の和に属する。
ここでは核上の固有値2が平行移動の存在を保証しているため、一般位置だけへの主張ではない。

## 6. rank4 接触配置：実一変数族と未被覆部分

Q=w+uv としてよい。\(\Omega=du\wedge dz+dv\wedge dw\)、
\(\lambda=2c-1\)、\(\mu=1-c\)、\(\rho=2-2c\) と置く。

\[
 Z=\lambda u\partial_u+\mu v\partial_v+\rho z\partial_z+cw\partial_w
\]

は \(\mathcal L_Z\Omega=\Omega\)、ZQ=cQを満たす。
(2) の complementary 係数を取った二次 lift は h を Casimir として

\[
\begin{aligned}
 P_{uv}&=\lambda hu,&P_{uz}&=-\mu hw-cuv,&P_{uw}&=-\lambda u^2,\\
 P_{vz}&=-\mu v^2-\rho hz,&P_{vw}&=-hw-\rho uv,&P_{zw}&=cvw-\rho uz.
\end{aligned}
\]

全 c 上で affine Jacobi は恒等的に零、
\(\operatorname{Pf}(P)=(c-1)(hw+uv)^2\)。従って c≠1 では所要の実際の rank4 Poisson 族である。
`family` mode はこの全パラメータ恒等式を確認し、c=2における全350×150射影Jacobi微分も扱った。
非零 modular minor は標数零の階数下界だけに使用した。
Jacobi微分102、鉛直24、軌道+scale46、これに真のc微分を加えて47であり、
全 Jacobi 微分が全ての実方向を零に写す恒等式も厳密に確認した。

従って射影接空間次元は高々23、実族の像閉包次元は少なくとも22である。
この不一致から成分性、一般被約性、既知成分への所属を結論しない。
証明書は `reports/CLASSIFICATION_TANGENT_DOUBLE_QUADRIC_FAMILY_20261107.json`、実行は約0.30秒。
共通 checker の `rank3` mode は別の rank3 調査用であり、本節の被覆証明に混ぜていない。

## 7. rank4 接触配置の全点を二つの源へ還元する

### 7.1 rank2 の Ω は存在しない

\(R=\ker\operatorname{Hess}Q\) は二次元、Q の一次部分 ℓ は R 上非零である。
Ω が rank2 のとき、第4節の X は固定二平面 K=kerΩ に値を持ち、
\(XQ=-kQ\)、\(\operatorname{tr}L=0\) を満たす（L=DX）。
Hessian の係数比較から R は L 不変で、商上の trace は −k。
従って R 上の trace は k。一方、一次係数を R に制限すると
\(\ell L|_R=-k\ell\) なので R 上の固有値は −k,2k である。
両方とも非零だから \(R=L(R)\subset\operatorname{im}L\subset K\)、従って K=R。
このとき商上の L は零となり、その trace −k≠0 と矛盾する。

### 7.2 symplectic Ω と Hessian の同時標準形

全点で Ω は symplectic。\(S=\Omega^{-1}H\) とすると
\([L,S]=(1-c)S\)、c≠1なので、全 \(\operatorname{tr}S^m=0\) から S は nilpotent。
S は symplectic Lie algebra に属し rank2である。
もし Jordan 型が (3,1) なら長さ3の鎖の端 S²a は、
\(\Omega(S^ix,y)=(-1)^i\Omega(x,S^iy)\) によって鎖と残る核ベクトルの全てに直交する。
非退化性に反するので S²=0、型は (2,2) である。
従って \(R=\ker S=\operatorname{im}S\) は Lagrangian 二平面。

symplectic 座標 x=(u,v), y=(z,w) と平行移動によって

\[
 \Omega=du\wedge dz+dv\wedge dw,\qquad Q=q_2(x)+b^ty,
\]

とできる。q2 は非退化二次式、b≠0。
L の block は \(\left(\begin{smallmatrix}A&0\\C&I-A^t\end{smallmatrix}\right)\)、C は対称。
二次項と y の一次項の比較は

\[
 A=\frac c2 I+T,\quad T\in\mathfrak{so}(q_2),\qquad Ab=(1-c)b.
\]

従って \((2-3c)q_2(b)=0\) である。

### 7.3 isotropic の全解は既約五変数源の像

q2(b)=0 のとき、Ω の標準形を保って Q=w+uv とできる。
全ての affine Z は、適当な p,q,r,e によって

\[
\begin{aligned}
 Z_u&=(2c-1)u+p,&Z_v&=(1-c)v+q,\\
 Z_z&=ru-qv+(2-2c)z+e,&Z_w&=-qu-pv+cw
\end{aligned} \tag{5}
\]

と書ける。逆に任意の c≠1,p,q,r,e に対し (4) は成立する。
従って全解は既約源 \((\mathbb A^1\setminus\{1\})\times\mathbb A^4\) の実多項式像で尽くされる。
一般の c、具体的には c∉{0,1/2,2/3,3/4,1} では、L の四固有値
2c−1,1−c,2−2c,c は相異なり非零である。
Z の唯一の零点へ移し、固有 symplectic 座標を取れば、第6節の一変数族へ移る。
Ω が対を作る固有空間の重みの和は1、Q の重みは c であることにより、
Q の二次部分と一次部分もそれぞれ uv と w へ規格化できる。
このとき二次部分は Hessian の radical R を通じて商へ降りている。
従って特殊な c で重みだけなら可能となる R 変数の二次項を新たな自由度に数えない。

従って (5) の全像は、一変数族を PGL で飽和した像の閉包に含まれる。
例外 c を個別に一般点と同一視せず、既約源の稠密開と閉性から全点を戻している。
独立証明書 `reports/CLASSIFICATION_TANGENT_DOUBLE_QUADRIC_T18_REVIEW_20261107.json` は、
c=2 の点 P2 に対する真の Poisson pencil P2+tR を構成し、P2+8R を T18 の実源
(p0,p1,t0,t1)=(10,12,−2,−2) と全25行列成分で同定している。
period と補正項の分母は共に8で非零である。
重み (0,1,−1,−2,0) の軌道極限が P2 を与えるので P2 は T18 に属する。
全350×150接空間行列の独立102次 minor による次元上界23と dim T18=23 が一致し、
P2 は全 Poisson 空間の滑らかな点で、その唯一の成分は T18 である。
従って既約 c 族の非空開は T18 に入り、その全体が閉性で T18 に入る。
上の既約五変数源もその閉包に含まれるので、isotropic 枝の全支持が T18 で被覆される。

### 7.4 anisotropic の全解は一つの軌道

q2(b)≠0 なら c=2/3。二次元直交 Lie algebra の元 T が非零ベクトル b を殺すので T=0。
従って A=I/3、I−Aᵗ=2I/3である。
L は可逆で、二つの異なる固有値1/3,2/3を持つ半単純な行列へ symplectic shear で移る。
Z の唯一の零点へ平行移動した後、重みの比較によって Q は x の非退化二次式と y の一次式だけとなる。
anisotropic の条件と線形座標・Q および Ω の尺度から、全て

\[
 Q=w+u^2+v^2,\qquad
 Z=(u\partial_u+v\partial_v+2z\partial_z+2w\partial_w)/3
\]

へ移る。従って残る点集合は一つの projective PGL 軌道である。
分母3を払った逆形式 \(\alpha=3Q\Omega-dQ\wedge\iota_{3Z}\Omega\) の対応する lift は

\[
\begin{aligned}
 P_{uv}&=hu,& P_{uz}&=-hw-3u^2-v^2,&P_{uw}&=-2uv,\\
 P_{vz}&=-2hz-2uv,&P_{vw}&=-3hw-u^2-3v^2,&P_{zw}&=4uw-4vz.
\end{aligned}
\]

全 Jacobi は厳密に零、Pfaffian は \(-3(hw+u^2+v^2)^2\)。
新しい `anisotropic` mode と
`reports/CLASSIFICATION_ANISOTROPIC_DOUBLE_QUADRIC_20261107.json` はこの点の全150方向を検査した。
Jacobi 階数下界99、鉛直24、軌道+scale47により、射影接次元は高々26、軌道次元は少なくとも22。
これは当該軌道の成分性を示さない。実行は約0.39秒。
この接空間の差を局所環の高次計算で埋める必要はなく、次の実変形でこの一軌道を被覆できる。

### 7.5 追加一軌道も MPS40 に属する

前節の P に、次の斉次二次 bivector R を加える。表示していない係数は零である。

\[
\begin{aligned}
 R_{hu}&=-hu/2,&R_{hv}&=-3hv/2,&R_{hz}&=hz/2,&R_{hw}&=3v^2/2,\\
 R_{uv}&=-uv/2,&R_{uz}&=uz/4,&R_{vz}&=vz.
\end{aligned}
\]

独立な有理数 sparse 外積計算により、全 t で
\(E\wedge[P+tR,P+tR]=0\) が成立し、その実退化五次式は

\[
 -\frac38(tuz+4hw+4u^2+4v^2)
 (thvw+tv^3+2h^2w+2hu^2+2hv^2).
\]

係数モノミアルの重みから bivector の二つの添字の重みを引く規約を使う。
(h,u,v,z,w) の重み (1,0,0,−3,−2) では P+R の全項の重みが非負で、
重み零の極限は R から \(R_{hw}=3v^2/2\) の項だけを除いた logcanonical bivector である。
その五次式は \(-3huvzw/8\) なので、実際に五独立超平面の SNC log symplectic 点である。

P+R を含む全 Poisson 空間の既約成分 C を一つ選ぶ。
連結な PGL5 は各既約成分を保つから、上の軌道極限も C に属する。
MPS Theorem1.4「SNC log symplectic locus と交わる全成分は既知40個」を適用して、C はその一つである。
さらに重み (0,1,1,2,2) に対し P の重みは −1、R の重みは0。
射影スカラーを含めると P+R の別の軌道極限が P を与えるので、P も C に属する。

この議論は単に「既知成分に属する特殊化があるから一般点も属する」と推論していない。
最初から P+R を含む同じ成分を選び、その PGL 不変性で SNC 点まで運んだ後、
SNC に交わる成分の全てを尽くす定理を適用する。
T33 の smoothing diagram と対応することも検査済みだが、全40成分の和への包含には
個別の T33 成分への唯一の所属は不要である。

証明書 `reports/CLASSIFICATION_ANISOTROPIC_MPS_DEFORMATION_REVIEW_20261107.json` は
SymPy と本ノートの producer を使わない Fraction 演算で全 t の Jacobi、五次式、両軌道の全重みを
約0.021秒で検査している。詳しい独立証明と正規形の査読は
`TANGENT_DOUBLE_QUADRIC_MPS_COVERAGE.md` に記録する。
これで第7節の全接触配置が既知成分和に被覆された。

## 8. 残る義務と監査境界

rank5 の全点認識、rank4 の空配置、rank4 非接触配置の全 β=0 帰着、
rank4 接触配置の全 MPS40 被覆が今回の支持結論である。
従って実退化因子 H+2Q、rank Q≥4 の全層が被覆された。
rank3 以下の inverse-form の追加クラス、他の退化五次因子型は本ノートの範囲外である。
特に「Q45 成分の閉包上にある」ことだけから全 affine ファイバーの被覆を推論していない。
この限定された幾何的被覆をもって、全 \(\mathbb P^4\) Poisson 成分数の上界が確定したとはしない。

## 9. 独立幾何査読

統合担当が旧smooth-double-quadric証明第3–7節を読み、新しいrank4の二配置を
上の行列計算とは独立に点検した。rank4頂点回避配置のrankΩ=2分岐では、
div X=0が定数二形式dβの核方向の成分から従い、非退化Hessianのtraceがk=0を強制する。
rankΩ=4分岐のtrace条件もc=1を強制し、全階数を尽くしている。
非接触配置では臨界直線上のZの固有値が正確に2となるため、一般点への限定なく
平行移動で一次項を消せる。Pの係数はαのcomplementary係数の定数倍なので、
αの定数+二次という性質はPにも移る。beta零への接続に欠ける例外点はない。

線形checkerの全90変数・140出力・20原始係数とQQ上の合成零・階数判定も読解した。
旧全150係数の成分証明をこの査読で再実行したとはしない。
第6節の族は実際の族として受理し、22と23の次元差を成分性の証明としない。
rank3の横断配置には別の全支持命題が得られ、
`RANK_THREE_QUADRIC_SUPPORT.md`とその独立reviewに接続する。
