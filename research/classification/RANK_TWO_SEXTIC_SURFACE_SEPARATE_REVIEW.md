# 次数六・rank 二支持の全 G 非優越：独立幾何・Fraction 監査

2026-11-07。reviewer は `log_manuscript`。
smooth-normalization producer `generic_exhaustion` の
`RANK_TWO_SEXTIC_SURFACE_NONCONICAL_DIVISOR_EXCLUSION.md` と、
conical-normalization producer `plane_collision_review` の
`RANK_TWO_SEXTIC_CONICAL_DIVISOR_NONDOMINANCE.md` を別に読解した。
指定された **一つだけ** の S(2,3) 射影を別の Fraction 実装で再構成した。
以下の全幾何監査と合わせ、両命題を受理する。

**受理する結論。** 実退化因子が既約被約 quintic である rank 四 Poisson tensor について、
次数六 integral 非正規曲面成分上の一般 restriction rank が二の全分岐は、
smooth normalization 二型では空、singular normalization 014/005 型でも空、
残る 023 型では全 PGL 飽和射影像が24次元以下である。
従って必要次元25--27の未知 G 優越成分の一般元にはならない。
E 成分帰属、全点被覆、他次数又は全 G/E の完了は含めない。

## 1. 全入力と計算の射程

* 六 file packet の入口：
  `inputs/scroll23_conductor_packet_manifest_20261107.json`。
  manifest 内五 file と manifest 自身を読み、全 hash を固定した。
* 新しい独立実装：`code/classification_rank_two_sextic_review.py`。
  producer、SymPy、CAS は import しない。
  reviewer 自身の既存 `classification_projected_veronese_zero_review.py` の
  標準ライブラリ Fraction・疎多項式・行列消去だけを再使用した。
* 契約：`jobs/rank_two_sextic_review.json`。
  新しい射影探索、点列、GB、高次 jet は行わない。
* 独立結果：`reports/CLASSIFICATION_RANK_TWO_SEXTIC_REVIEW_20261107.json`。
  再 gate 後0.083秒、最大RSS13,952 KiB、全項目 PASS。

初回は、implicit substitution の非零 support 行数を全 \(H^0(5H)\) の81 basis 行数と
同一視した checker の寸法 assert で、report 生成前に停止した。
全81 basis 行（この処方で零となる行も含む）を明示して修正した。
producer bytes、射影、幾何の仮定又は成功条件は変更していない。
失敗した実行は証拠として使わない。

## 2. rank に依存しない正規化と全 node/cusp

既受理の length/Bertini/Höring の部分だけを使用する。
\(\deg S=6\) は \(g=0,\delta_S=1\) を強制し、二次元の非正規支持は S 一つ。
\(X\) の完全 \(|L|\) 埋込みは P7 の minimal-degree threefold であり、
全型は 113、122、014、023、005 である。
rank 零用 cone 排除は今回の rank 二に流用しない。

まず 005 は他の型の class-group 計算を使う前に処理出来る。
line vertex は有限線形射影で同次元の vertex に写る。
像 D は plane quintic curve とその line の join である。
vertex 外は base curve と affine space の積だから、二次元非正規支持は
base curve の singular point と line の join、すなわち plane の有限和である。
vertex 自体は一次元であり、次数六 integral S は存在しない。

113、122 は smooth projective bundle。
014、023 では \(\widetilde X=\mathbb P(O\oplus O(a)\oplus O(b))\to X\) が
zero-summand section P1 だけを点へ縮める small resolution である。
exceptional locus は codimension 二、divisor は一つも縮めない。
したがって \(\operatorname{Cl}X\simeq\operatorname{Pic}\widetilde X=\mathbb ZL\oplus\mathbb ZF\)、
divisorial conductor の strict transform は
\(C\sim-K_{\widetilde X}=3L-3F\) である。
vertex の導手 ideal が Cartier 又は CM だと仮定する必要はない。

\(\delta=1\) の全 generic transverse type は node/cusp。
cusp は \(C=2\Gamma\) となり、上の Pic の L 係数三に反する。
よって残る四型では generic node、C は reduced であり、
一つの prime（二重 cover、非分裂 node を含む）又は二つの prime の全場合を持つ。

## 3. node の留数非零性と fixed D の全 tensor 一意性

完成 node 環 \(K[[u,v]]/(uv)\) の各 minimal prime は任意の Hamiltonian derivation に保存される。
各 normalization branch は全環の surjective 商なので bracket は regular で、branch の conductor ideal も保存する。
係数体 K は Poisson 中心と仮定しない。
既受理の completion continuity・有限分離拡大・faithful descent が同じく適用出来る。
よって \(\pi_X\) は全 codimension 一で regular かつ conductor に tangent。

inverse symplectic form の closed residue を \(\alpha=i_{\pi_X}\Omega_X\) とする。
conductor 一般点の \(\Omega_X=(du/u)\wedge ds\wedge dt\) で、
tangency により mixed normal bracket は u で割れ、
\(\alpha=c\,du/u+\text{regular}\) となる。
restriction rank 二と conductor から S への generic separability により \(c\ne0\)。
closedness は各 prime 上の residue を非零 constant にする。

single prime なら \(L^2\Gamma=12\)。一般二 hyperplane の smooth pullback curve は
SingX、vertex、全 codimension 二 bad locus を避ける。
residue theorem が \(12\lambda=0\) を与え、rank 二からの \(\lambda\ne0\) と矛盾する。
非分裂 node の二局所 branch が一つの大域 prime を作る場合も、この一つの constant で排除する。
smooth 型では residue 零から \(H^0(\Omega_X^1)=0\) を使う producer の論証も有効だが、
上の非零性だけでも同じ矛盾を得る。

split なら \(C=\Gamma_1+\Gamma_2\)、各 \(\Gamma_i\to S\) は degree 一なので \(L^2\Gamma_i=6\)。
全 residue pair は \(6\lambda_1+6\lambda_2=0\) に入り一次元以下である。
同じ residue の closed forms の差は全 codimension 一で holomorphic。
X の smooth locus へ、さらに small resolution の exceptional curve 外へ引き戻し、
smooth \(\widetilde X\) の locally free \(\Omega^1\) の Hartogs で全体へ延ばす。
rational P2-bundle over P1 では \(H^0(\Omega^1)=0\) なので、この差は零。
従って fixed D の全 \(\alpha\) 空間は一次元以下である。

dualizing generator \(\Omega_D\) を固定した contraction は generic smooth locus 上で bivector に injective。
ambient restriction も \(H^0(\mathbb P^4,\Lambda^2T(-5))=0\) で injective。
actual Pfaffian の scalar が各 \(\pi\) で変わっても contraction form は closed residue の constant 倍なので、
この一意性は同じ D を持つ全 projective tensor に適用する。
すなわち fixed-D projective fiber は高々一点である。

## 4. smooth 二型と 014 の split conductor を尽くす

各 prime \(\Gamma_i\sim a_iL+b_iF\) は一般 P2 fiber 上で effective なので \(a_i\ge0\)。
\(a_i=0\) の prime は fiber F だけだが次数一なので degree 六に反する。
よって \(a_i\ge1\)、和は三、各次数は六から、交換を除く唯一の可能性は
\[
 \Gamma_1\sim L+F,\qquad\Gamma_2\sim2L-4F.
\]
全 \(|2L-4F|\) section を fiber 座標 u,v,w で書くと、

* 113：\(w(au+bv+q_2(s,t)w)\)。
* 122：\(av^2+bvw+cw^2\)。
* 014：\(w(au+\ell_1(s,t)v+q_4(s,t)w)\)。

最初と最後は非零 fixed divisor w と他の非零 divisor の積、122 は C 上で binary quadric が二一次因子へ分解する。
double factor も prime effective divisor ではない。
small resolution は divisor を縮めないので、014 の因子を vertex へ消す余地も無い。
single-prime は §3 で既に除いたため、三型全てが空となる。

## 5. 一つの指定 023 射影の全 conductor を独立再構成

\(Y=S(2,3)=\mathbb F_1\)、\(H^2=5,Hf=1,K_Y=-2H+3f\) とする。
指定四 section の係数行は
\[
 A=(s^2,st,t^2,0),\qquad B=(s^3,t^3,s^2t,s^3+t^3).
\]
全六 minor を別に展開し、s一での univariate gcd 一、s零での二行 rank 二を確認した。
従って複素 P1 の全域で rank 二、basepoint は無い。
全7行4列 coefficient matrix の rank 四も再構成した。
H ample より射は finite、四 section 独立より像は plane でない。
\(5=\deg\nu\deg T\) の素数性から \(\deg T=5,\deg\nu=1\)、Y は像の全正規化である。

新実装は resultant や producer 五次式を入力の解として使わず、
全 \(H^0(5H)\) の81行と四変数次数五の56 monomial から代入行列を構成した。
rank は55、独立55minorは非零で、primitive な一つの kernel relation が保存 quintic と完全一致した。
homogeneous pullback の全係数零も確認した。
像の integral 性と次数五から、この式の絶対既約性も独立に従う。

s=v一 chart で \(d=1+T^3\)、
\(Z=(r+1,rT+T^3,rT^2+T,d)\) とする。
quotient derivative を別に計算した Jacobian numerator は
\(J_{\rm num}=-T^3r-3T^3-3T^2-r\)。
再構成五次の homogeneous z derivative を全て代入し、多項式除算で
\[
 F_{Z_2,\rm hom}(Z)=J_{\rm num}G
\]
を余り零、scalar 一で得た。
これは \(G=dF_z/J\) であり、必要な d 因子を省略していない。
G の全係数は保存値と一致した。

有限双対性から \(cO_Y=\omega_Y\otimes\nu^*\omega_T^{-1}\) は可逆層で、
その実導手 divisor は \(3H-3f=3C_0+6f\) である。
canonical frame \(dx\wedge dy/F_z\) は \(O_T(1)\) の Z3 section に対応する。
\(T'=1/T,r'=r/T\) の \(dT'\wedge dr'=-T^{-3}dT\wedge dr\) と、
\(\operatorname{div}(d)=F_d-3F_\infty\)、\(\nu^*Z_3=C_0+F_d\) を使うと
\(C=\operatorname{div}(G)+3C_0+6F_\infty\) となる。
全項 \(r^iT^j\) の正しい global monomial は \(u^iv^{3-i}s^{6-i-j}t^j\)。
独立展開は global section と一致し、境界値も
\[
 G|_{s=0}=-t^4u^2v\ne0,\qquad
 G|_{v=0}=-su^3(4s^2+6st+3t^2)\ne0
\]
となった。d零も polynomial identity により復元され、missing component は無い。

squarefreeness は別の短い証明で確認した。
G を \(\mathbb Q[T][r]\) とした coefficient content は一。
T零は degree 三を保ち、\(G(r,0)=-4r^3-6r^2-3r-1\) とその derivative の gcd が一。
したがって \(\mathbb Q(T)[r]\)、さらに primitive な \(\mathbb Q[T,r]\) で squarefree である。
char 零なので C 上でも reduced となる。

全31未知数の closed-form 条件
\(G(A_T-B_r)-AG_T+BG_r=0\) を別に全44行として構成した。
全 matrix entry、row/column labels、\((G_r,G_T)\) kernel が保存値と一致する。
30次 minor は \(3486368987873280\ne0\)、よって C 上でも kernel は一次元。
各絶対因子から \((G/G_i)dG_i\) が全次数制限内に入り、その \(d\log G_i\) は prime residue で独立。
従って絶対因子は一つしかない。
境界二 divisor の非零 restriction と合わせ、**全 global conductor が integral** と従う。
一つの Poisson tensor の標本でも、QQ 因数分解を C へ流用した論法でもない。

## 6. 全有限射影における相対導手の base change と properness

023 の全五 section は basepoint-free open \(U\subset\operatorname{Gr}(5,8)\)、次元15。
vertex で非零なので \(W=V\cap H^0(Y,H)\) は四次元。
W に basepoint があればその coneline 上で最後の section も零に出来て V の basepoint となる。
逆向きも成立するので U は basepoint-free
\(U_4\subset\operatorname{Gr}(4,7)\)（irreducible、次元12）上の三次元 affine bundle である。
fiber は \(H^0(Y,H)/W\)。vertex translation により D は四 section の surface image の cone となる。
全 W の像は上と同じ素数次数の議論で finite birational quintic である。

universal degree 五 restriction map の kernel は全 W で一、rank 一定なので line bundle。
その式が universal P3-bundle 内の flat relative quintic hypersurface \(\mathcal T\) を作る。
\(Y\times U_4\to\mathcal T\) は proper かつ全 point fiber finite なので finite。
全 fiber が正規化であるため、これは total space 上でも finite birational。
source は smooth、target は相対 hypersurface Gorenstein である。

**導手が fiber ごとの実導手であることも確認する。**
局所の target A は Gorenstein、有限 source B は smooth かつ同次元なので maximal CM A-module。
従って \(\operatorname{Ext}_A^{>0}(B,A)=0\)。
base の regular parameter t に対して A と B は flat、t は両者で非零因子である。
\(0\to A\xrightarrow{t}A\to A/t\to0\) へ Hom(B,−) を適用すると
\[
 \operatorname{Hom}_A(B,A)/t
   \simeq\operatorname{Hom}_{A/t}(B/t,A/t)
\]
となり、全 base parameter を順次用いて fiberwise finite duality と一致する。
evaluation at 1 による Hom の B への inclusion も自然で、この base change と可換である。
よって候補の limiting cycle ではなく、各 fiber の conductor ideal 自体を与える。

relative duality は source 上で
\(\omega_{Y/U_4}\otimes\nu^*\omega_{\mathcal T/U_4}^{-1}
 =O_Y(-3H+3f)\otimes\text{base line}\)
となる。inclusion は各 fiber で非零であり、base line を局所 trivialize して
\(U_4\to\mathbb P H^0(Y,O_Y(3H-3f))\) という morphism を与える。
normalization defect の任意の flatness を無証明で課したのではない。

F1 の effective cone は \(C_0,f\) なので、\(3C_0+6f\) の二つの非零 effective class への分解は
\(0\le a\le3,0\le b\le6\) の有限個。
各 section multiplication の projective image は proper で閉であり、その有限和が全 nonintegral divisor locus。
§5 の一 witness が integral なので、その inverse image は \(U_4\) の proper closed subset、次元11以下。
split conductor locus はこれに含まれ、U へ戻して14以下となる。
vertex 外の polynomial extension は normalization/conductor と可換だから、
surface conductor の splitting と X の divisorial conductor の splitting は同じである。
vertex 点から divisor が追加されることは無い。

## 7. 全像の次元と最終判定

023 の実 polarized automorphism subgroup は
\(E=O\oplus O(2)\oplus O(3)\) の \(h^0(\operatorname{End}E)=12\)、
scalar 一を除き base PGL2三を加えた14次元。
small resolution の dense open で faithful なので X でも faithful である。
target projective frame は24、全 split framed-map source は高々 \(14+24=38\)。
precomposition は D を変えず、finite birational map を固定する automorphism は function field 上で identity、
従って全体で identity。よって source fiber は実際に14次元で、D 全像は24以下となる。
§3 の fixed-D projective Poisson fiber は高々一点なので、全 Poisson 像も24以下。
ここでは morphism \([\pi]\mapsto[\operatorname{Pf}\pi]\) の rank 四 open 上の algebraicity と fiber-dimension 不等式を用いる。

以上で smooth113/122、cone014/005 の空性、cone023 の全像24以下を別々に確認した。
全 node/cusp、single/non-split/split、全射影と特殊 codimension 二導手を保持する。
両 producer 命題を合わせた **全 degree 六 rank 二 G 非優越**として受理する。
degree 五以下、plane の他の zero divisor、全 E 帰属、全分類又は新成分数を閉じたとはしない。
103成分・152頁の保存済み主TeX/PDF、global 台帳、Git は今回変更しない。
