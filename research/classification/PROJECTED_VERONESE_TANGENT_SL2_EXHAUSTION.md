# projected Veronese全十五源：SL2 stabilizerからactual既約quinticを排除する

2026-11-08、PRIMARY E-FINITE-SURFACE。
**全tangent Poisson sourceのrank四・actual被約既約quintic枝は空。**
sourceの選択点、Jacobiの追加operator、GB、minor、jet又はprofile列を使わない。
既存の全source次元上界と、一回の完全SL2 module証明書を結ぶ。
新componentの増分は零。別作者によるmodule再生は別receiptで管理する。

## 1. 使う二つの既受理全源定理

S⊂P4をsmooth projected Veronese surfaceとする。
全て一つのPGL5軌道なので、

\[
 q(a,b,c)=(a^2-b^2,\ b^2-c^2,\ ab,\ ac,\ bc)
\]

へ正規化してよい。
`PROJECTED_VERONESE_TANGENT_NONDOMINANCE.md` と同独立reviewにより、
全global bivectorをSへtangentにするEuler商の空間Vは正確に十五次元。
canonical divergence-zero liftで一意に表され、全350 Jacobiはhomogeneous quadratic equationsである。
その**全Poisson cone** X⊂Vのaffine次元は二以下。
この上界は全54生成元・242polynomial・10,079exact liftsの受理済み証明を使う。
今回再びGBを走らせたものでも、chosen slice又は一つのtangent rankからの上界でもない。

また `PROJECTED_VERONESE_ZERO_SURFACE_REVIEW.md` により、S上で零となる全canonical sourceは五次元で、
その中のPoisson tensorは零だけである。
以下ではこの**全rank零source定理**を明示入力として用いる。

## 2. 実際のSL2作用と完全module

非退化二次式a²+b²+c²を保つSO3がprojection centerを保ち、S及びambientP4に線形作用する。
u=a+ib,v=a−ibとすれば、SL2のSym²表示の三fieldを

\[
 H=2u\partial_u-2v\partial_v,\quad
 E=u\partial_c-2c\partial_v,\quad
 F=2c\partial_u-v\partial_c
\]

と取れる。SL2の有限中心を除いた実image G≅PGL2は三次元の連結代数群である。
これは単に形式Lie代数作用を仮定したのではなく、projectionのcenterを保つ実線形群作用。
全canonical lift、tangency、Jacobiはこの群で保たれる。

ambient座標x0,…,x4の次の実可逆線形変更を使う。

\[
 y_0=x_0+2ix_2,quad y_1=x_3+ix_4,quad y_2=x_0+2x_1,quad
 y_3=x_3-ix_4,quad y_4=x_0-2ix_2.
\]

これらはq上でu²,uc,uv−2c²,vc,v²となり、Hのphysical coordinate weightsは

\[
                             (4,2,0,-2,-4).             \tag{1}
\]

全十五basisに対しpushforward infinitesimal action
ρ(A)P=AP+PAᵗ−(Ax)·∇Pをそのまま展開した。
全150×15係数の三actionを確認し、全chain matrixと両側inverseから

\[
                    V=\operatorname{Sym}^2\oplus
                      \operatorname{Sym}^4\oplus
                      \operatorname{Sym}^6              \tag{2}
\]

を厳密に得た。全てSL2標準二次元表示の対称冪であり、multiplicityは一。
Hのcharacteristic polynomialだけから(2)を推測したのではない。
最低weight−2,−4,−6の各E-kernel lineからFで長さ3,5,7の全chainを作り、
合計十五列の行列とinverseの両積を確認した。
source conventionは[H,E]=−2E,[H,F]=2F,[E,F]=−Hで、Eがlowest-weight lineを消す。

さらに、全rank零五basisは(2)のSym4 blockに正確に等しい。
全三fundamentalfieldのwedgeのEuler商はSym2 blockに正確に等しい。
それぞれ5×5、3×3のblock inverseと他blockの全零を保存し、全sourceの等号を確認した。
三つの選んだtensorを全supportと見なしたのではない。

## 3. 全Poisson点がBorel eigenlineへ入る理由

Jacobiはhomogeneousなので、非零Poisson tensorの射影支持P(X)はdimension一以下。
Gは連結なので有限個の既約componentをそれぞれ保つ。
従って任意の点[Π]∈P(X)のG-orbitはdimension一以下で、

\[
                    \dim\operatorname{Stab}_G[\Pi]\ge2.
\]

PGL2のdimension二の連結真部分群はBorelであり、dimension三なら全群。
いずれも適当な共役Borelを含む。
従ってあるg∈GでgΠのlineが標準Borelに保たれる。
unipotent radicalはGmへの非自明な代数的characterを持たないので、
そのgenerator EはgΠを零にし、torus Hは同lineにscalarで作用する。

(2)の最低weightは相異なるため、全Borel eigenlineはSym2、Sym4、Sym6の
三lowest-weight lineのいずれかだけである。
同じweightのsummandをmixする自由parameterは無い。
この全支持の矢印はprior **全cone dim≤2**から得たものであり、
normalformの点でJacobiを確かめて全supportを外挿したものではない。
原点はrank四でなく、schemeの非被約厚みは幾何閉点を追加しない。

## 4. 三つの強制branchはrank四のactual既約quinticを持たない

**Sym2。** 全spaceはfundamentalfieldのwedge像である。
三次元Lie algebraのΛ²の任意の元は二vectorのwedgeなので、
対応するprojective bivectorのrankは二以下。
canonical liftにEuler補正があっても射影rankは変わらない。
よってこの全orbitにrank四点は無い。

**Sym4。** 全spaceがS上rank零の五源に等しい。
先行の全rank零PV theoremにより、Poissonなら零であり、非零射影点は無い。
rank二restrictionに対しrank零源を無断で仮定したのではなく、完全moduleのこのbranchだけに適用した。

**Sym6。** 最低weight representativeΠはρ(H)Π=−6Πを満たす。
SL2 ambient表示のdeterminantは一であり、actual quintic Fを
E_ambient∧Π²=F·volumeとして定めると、Pfaffian mapはequivariantな二次mapである。
torus tの作用から

\[
 F(t^{-4}y_0,t^{-2}y_1,y_2,t^2y_3,t^4y_4)=t^{-12}F(y).
\]

従ってFの各degree五monomialは(1)についてweight十二を持つ。
y0を含まないdegree五monomialのweightは高々2·5=10なので、

\[
                              F=y_0F_4.                \tag{3}
\]

F=0ならrank二以下。非零ならactualDはproper hyperplane factorを持ち、被約既約quinticに反する。
反対のBorel conventionならweight−十二となりy4因子を得るだけである。
全G共役と非零tensor scaleでも因子性は保たれるため、この全orbitが排除される。
ここでactualPfaffianを新しい点へ代入する必要は無い。

以上で全三branchを尽くし、全tangent十五源のrank四・actual既約被約quintic空性が従う。
Sが非正規曲面componentである追加仮定すら、この強いtangent-source排除には不要。

## 5. 分類へのscopeとexact certificate

全smooth projected Veroneseは一つのPGL5軌道なので、全埋込みへ結論を戻せる。
actualDが可約なrank四Sym6点の存在／membership、全E又は総component数の完成は主張しない。
それらは既存の最上位factor義務へ戻す。新componentの増分は零。

有限契約 `jobs/projected_veronese_tangent_sl2_exhaustion.json` のgate後、
**4.374286秒／RSS63,560KiB、全module identity PASS**。
全150×15 action、座標left inverse、全十五chainとinverse、zero5/wedge3 block inverse、
physical weight座標とinverse、embeddingとcenterの全恒等式、入力SHA、inline sourceは
`CLASSIFICATION_PROJECTED_VERONESE_TANGENT_SL2_EXHAUSTION_20261108.json`へ固定した。
旧dimension証明のGB又はJacobi operatorを再実行していない。
既存source、oldreport、currentjob、主稿、global台帳及びGitは変更していない。
