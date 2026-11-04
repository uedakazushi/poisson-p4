# v=0からの小さいactual β sectionとfull tensor lift

2026-11-04。元の六sourceのv=0境界から、同じactual B_NEW15へdenseな新しい有理sectionを
作る。第一積分の約分・逆変換だけでdominanceを推測せず、元の多項式kernelの横断性と
DVRのvaluationを用いて実際のrawjet像へ結び付ける。

新sourceは**非零β座標9個、分子の総項数12、最大3項、全分母monomial**となる。
旧t=2 sourceのh=s²−25rq²は分母に現れない。旧sourceと新sourceの同じ見かけの座標名を
同一視しないため、新しい五変数をxi0,r0,s0,q0,B0と記す。
高次不変式・source26・新しい点やrankのjobは行わない。

## 1. 五つの定数をv=0へ厳密に制限する

A=256、alpha=384X+(3072W−1)、beta=384X−(3072W−1)とし、flow座標を
rf=A+v alpha、sf=A+v betaとする。これは後のshape r0,s0と違う座標である。
既存の第一積分で見かけの0/0を先に約分すると、全有理式の恒等式として

\[
U=u,\quad
K=\frac{\alpha(2A+v\alpha)s_f^2}{\beta(2A+v\beta)r_f^2},
\]
\[
I_v=\frac{r_fs_f-8uv(\alpha-\beta)^2}{s_f^2\alpha(2A+v\alpha)},\qquad
I_c=\frac{c r_f^3}{s_f\alpha^2(2A+v\alpha)^2},
\]
\[
I_b=\frac{B r_f^4}{s_f^2\alpha^3(2A+v\alpha)^3}.
\]

sourceはこの全五式と保存済みの元第一積分との一致を検査する。v=0のgeneric点では
rf=sf=A≠0、alpha beta≠0なので全てregularで、制限は

\[
\boxed{U=u,\quad K=\alpha/\beta,\quad
I_v=1/(2A\alpha),\quad I_c=c/(4\alpha^2),\quad I_b=B/(8A\alpha^3).}
\]

逆は有限な式

\[
\alpha=1/(2AI_v),\quad\beta=\alpha/K,\quad
u=U,\quad c=4I_c\alpha^2,\quad B=8AI_b\alpha^3,
\]
\[
X=(\alpha+\beta)/768,\qquad
W=(1+(\alpha-\beta)/2)/3072
\]

であり、両方向の全五合成を厳密に確認する。
conic parameterも先に約分するとt=A(alpha−beta)/(sf alpha)で、v=0ではt0=1−1/K。
このboundary sectionを旧t=2 sectionと同じ点だとは呼ばない。

tauはu,c,Xを反転しW,Bを保つため、(alpha,beta)→(−beta,−alpha)。上の制限は
既存の(U,K,Iv,Ic,Ib)上のtau作用と完全に可換である。
従って五定数体F5からv=0 source体k(u,c,X,W,B)への制限は**tau-equivariant birational**。
ただしこの事実だけではactual β像のdominanceはまだ従わない。

## 2. 例外的なblowdownを除くDVR論証

使うvector fieldはδ=(v/3072)Vではなく、元の保存済みの**多項式kernel V**である。
実際に全26 kernel成分がu,v,X,W,Bについて多項式で、perspective4成分は零、さらに

\[
\boxed{V(v)|_{v=0}=402653184\ne0}
\]

を確認した。c成分はdc=c V(log c)で、c≠0上regular。Darboux frameの分母はcBだけで
vを含まない。uもsource chartで単元である。

GL4×source上のraw value＋linear jet写像をPhiとし、そのimage closureをYrawとする。
保存されたraw30の全kernel恒等式を使い、source側Vとgdot=g A(source)というGL4側の
補償からderivation Dを作る。DはPhiの全30座標のpullbackを殺す。
v=0のgeneric DVRで、DはregularかつDvはunitである。

もしv=0 sourceのGL4 saturationがYrawの真の閉部分Zへ入れば、非零f∈I(Z)を取れる。
PhiはYrawへdominantなのでfのpullbackは非零であり、このDVRでv^n h、n≥1、h unitと書ける。
しかし

\[
D(v^nh)=n v^{n-1}(Dv)h+v^nDh
\]

のvaluationは正確にn−1で、標数0のnとDv,hがunitなので零でない。Df=0に矛盾する。
よってv=0のraw GL4 saturationはdense。
値のDarboux正常化・perspective除去はgeneric v=0でも定義されるため、そのprimitive βの
**Sp4 saturationはactual B_NEW15にdense**である。新しいminorや点の外挿を使っていない。

## 3. 新しい五変数とactual quotientのindex1

v=0上で

\[
r_0=u^2/1600,\quad s_0=48uX,\quad
\xi_0=40c/(Bu),\quad q_0=3072W-1,\quad B_0=B
\]

と置く。逆はu²=1600r0、X=s0/(48u)、c=xi0 B0 u/40、W=(q0+1)/3072。
従ってv=0 source体はF0(u)、F0=k(xi0,r0,s0,q0,B0)、tau u=−uであり固定体はF0。

既知actual β quotientはF5^tauと同じpullbackでindex1である。§1のtau可換なbirational
制限に加え、§2のactual β dominanceがあるので、真のinvariant rational functionsの
pole locusにsource全体が落ちる問題はない。従ってactual quotient体はこの新F0と同定される。
第一積分を単に境界へ延長しただけの主張ではない。

旧t=2のxi,r,s,q,Bと新xi0,r0,s0,q0,B0は別の関数である。
r0だけは旧rと共にu²/1600なので一致する。
他の座標変換は「新seedから§1の五定数を計算し、保存済みのfull rational inverseをt=2で
評価して旧seedを得る」という有限recipeで完全に指定する。
逆には旧t=2 seedの五定数を約分式で計算して、§1のv=0 inverseを使う。
これらの元のinverse文字列と新旧sourceへの全代入規約をcertificateに保持する。
大きいcoordinate-change係数へ展開したものではない。

## 4. 実際の16β係数：九つ非零、十二項

元の保存済み全16βへv=0を代入し、既受理のfinite pullback
g_u=diag(u,u⁻¹,1,1)を適用する。元のCartan weight w1によりβiをu^w1_i倍し、
u²=1600r0のquadratic étale algebraで全16のeven/odd部を計算する。
**全odd部が厳密に零**となる。

元W16のzero-based coordinate順で、非零成分は

\[
\begin{array}{c|c}
i&\bar\beta_{0,i}\\\hline
0&-45s_0/(16B_0^2\xi_0)\\
1&-30/(B_0^2\xi_0)\\
2&2400r_0\\
5&(5q_0+s_0+1)/(240B_0)\\
6&s_0\xi_0/(54000r_0)\\
8&-32B_0r_0\xi_0/3\\
9&-32B_0^2\xi_0^2/675\\
14&(s_0+1)/(120B_0)\\
15&-64B_0r_0\xi_0/3
\end{array}
\]

で、他の七成分は零。分子項数の合計12、最大3、全分母はr0,xi0,B0と有理定数のmonomial。
beta map自体にh0=s0²−25r0q0²の反転は不要。
§1のbirational constant chartではalpha beta≠0、従ってh0≠0を使うが、それをβ式の分母と
取り違えない。β式はさらに広いr0 xi0 B0≠0上で定義される。

新sectionのSp4 saturationはdenseで、actual invariant-field pullbackはindex1、generic β
stabilizerもidentityなので、Sp4×Spec F0→B15は既存sectionと同じ論法でbirationalとなる。
新しいsource点やrank jobは不要である。

## 5. 実際のfull tensor section

certificateにはP(u,1,0)の全affine式、点p0、Darboux R0、alpha0を保持する。
u²=1600r0上でc=xi0 B0 u/40を掛け、full tensorへ次のphysical変換を行う。

\[
\operatorname{block}(g_u^{-1})\,
\operatorname{Persp}(\alpha_0)\,
\operatorname{block}(R_0)\,\operatorname{Translate}(-p_0),
\quad
\alpha_0=(0,0,(s_0+1)/(240B_0),-32\xi_0B_0r_0/3).
\]

βのpositive-weight pullback g_uは、full tensorではphysical g_u⁻¹に対応する。
source pointのZは元の−3/(4u)(B+W/1536−W²+X²/64+2/9)へ上記代入をしたもの。
既受理の全次数tau同値とR_tau D=HR、g_−u=Hg_uにより、この**full tensor全体**はtau固定。
Euler quotientとS61.normalizeを施し、quadratic étale algebraのeven部を取る有限recipeを得る。

全factorとそのinverseはuとxi0,r0,B0のLaurent ringにあり、

\[
\det R_0=9/(8r_0\xi_0^2B_0^3)
\]

はunit。したがって下降後の全係数は
k[xi0^±1,r0^±1,B0^±1,s0,q0]に入り、h0分母は生じない。
今回はfull150係数へ再展開せず、この全入力を持つ有限operatorとfull equivarianceを使う。

v=0 normalized boundary像をY15とする。Sp4×Spec F0からfull tensor liftを作ったmapを
B15へ射影すると、前節のbirational mapになる。よって同じ函数体の両端に挟まれ
Y15→B15もbirationalである。従ってこのliftはgeneric βに対する**rational full tensor section**。
別の数え上げでも、v0 parameter代表はu±の二つ、full stabilizerはGaのみ、y0切断は一回なので、
元v0 full sourceはdegree2で、そのtauを除けばdegree1となり整合する。

任意のgeneric βからの有限recipeも指定できる。既存の旧quotient逆recipeで旧E5座標を回収し、
§3のconstant mapで新F0へ移す。beta_bar0を作り、既存のC/D intertwiner＋odd β scaleによる
平方根なしの一意Sp4 inverseを使う。そのphysical inverseを上のfull liftに掛ければよい。
generic freenessとdense sectionにより同じopenで実行できる。

この境界sectionのbirational性は、一般Y16への旧source degree8、tau下降後degree4を変更しない。
boundaryの特殊fiberを一般fiberの点数へ外挿しない。

## 保存と範囲

`code/q_new_family_v_zero_beta_section_20261104_v63.py`と同名report directoryに、
全約分式・両逆・Vv0・16β・項数とunits・新旧field変換・fulltensor operatorを保存する。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_v_zero_beta_section_20261104_v63.py \
  --output work/vzero-beta-section-replay
```

上限60秒/2GB/32MiB。NEW定義式の発見、degree26以降のsource/invariant job、全成分尽くし、
全β点でのsection regularity、full generic-fiber座標のdegree1は主張しない。
