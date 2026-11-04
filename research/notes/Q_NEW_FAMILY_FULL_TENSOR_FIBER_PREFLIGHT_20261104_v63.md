# NEW の full S61 tensor：有理な一次元 generic fiber への有限構成

2026-11-04。対象は実際の既約族の**被約閉包** C16⊂S61と、そのprimitive firstjet像
B15である。有限operator recipeによってfull quadratic tensorを一変数の有理体へ
持ち上げ、同じB15の函数体上でLürothを適用する。βのkernelをfull tensorの対称性と
同一視しない。結論の概念的独立監査は別packetとし、全projective成分25の有理性への
帰結は本稿に含めない。

## 受理済みの入力

1. 元の全三parameter族P(n0,n1,n2)は全次数でPoissonで、その全canonical150係数がある。
   n1=1とmodular Gaで評価点y=0へ移す操作はfull tensorの有限同値である。
2. 六source F6=QQ(u,v,c,X,W,B)からvalueとlinear jetを取り、全source6＋GL4 16＋
   perspective4のraw30×26行列はgeneric rank25。唯一のkernelはsource成分非零で、
   **保存されたperspective補償は全四成分が恒等的に零**である。
3. そのsource derivation δは積分済みで、F6=F5(t)、kerδ=F5。
   F5=QQ(U,K,Iv,Ic,Ib)であり、逆変換は全て有理式である。
4. 既知tauはu→−u、t→t/(t−1)と作用する。t=2 sectionの固定体は
   E5=QQ(xi,r,s,q,B)、F5=E5(u)、u²=1600rである。ここでshapeのr,sはflow座標の
   r,sと異なる。元のsource U=uはδの定数である。
5. 実際のB15はSp4×Spec E5とbirationalで、generic geometric stabilizerはidentity。
   これは単なるLie stabilizer零ではなく、C,Dの全共通commutantがscalarである証明を使う。
6. 元source＋Sp4から全S61の61座標への有限有理mapはC16にdenseである。
   既存の像次元16、β像次元15の証拠は再実行しない。

依存する原証明・独立gateのpath/hashはreportに固定した。今回の計算は以下の小さい
matrix identityと座標恒等式だけで、全150係数の新規展開や追加点を用いない。

## raw jet の族を同じ有限GL4 orbitへ合わせる

full tensorの評価点をp=(X,0,Z,W)、

\[
Z=-\frac3{4u}(B+W/1536-W^2+X^2/64+2/9)
\]

とする。平行移動と既存の有理Darboux行列R_tで値を標準Jへ移し、perspectiveを
まだ除去しないraw linear jetをL_tと書く。trace covectorは

\[
\alpha_t=\left(0,0,\frac{48Xu+1}{240B},-\frac{4cu}{15}\right)^{\mathsf T}.
\]

特に最後の成分は既知source chartで非零である。R_tは全六sourceを保った既存式であり、
t=2の参照raw jetを(J,L2)、covectorをalpha2とする。

**βだけでなくpair(β,alpha)のorbitを使う。** raw kernelのGL4補償をAとすると
δC+A C+C Aᵀ=0である。値RCRᵀ=Jを微分すると、正常化後の補償

\[
\mathcal B=\delta R\,R^{-1}-RAR^{-1}
\]

はsp4に入る。従って正常化後のraw pair(β,alpha)も同じsp4方向に沿って動く。
全pair-invariantの引戻しはkerδ=F5に入り、さらに後述の**full** tau同値があるので
F5のtau固定体E5に入る。一方、primitive βのactual invariant fieldは既受理のE5であり、
pair-invariantの部分体である。よってこの実際のpair像のquotient fieldも正確にE5である。
βのquotientだけからalphaが合うと推測していない。

t=2でg_uの正weight β gaugeも掛けた参照pairのcovectorは

\[
\bar\alpha2=(0,0,(s+1)/(240B),-32\xi Br/3)^{\mathsf T}
\]

である。これはu²=1600r、X=s/(48u)、c=xi Bu/40を代入した新しい小さい恒等式でも
確認した。従ってt=2のpair sectionそのものがE5上で有理である。

**orbitの存在はODEの形式的積分から主張しない。** 保存されたraw kernelがperspective零
なので、sourceのt方向はGL4のorbit方向に接する。raw source像のGL4 saturationをYとする。
GL4＋source部分のgeneric rankは21である。これは全26列の唯一kernelがperspective成分零
であるためで、GL4 actionのgeneric orbitは16次元、rational quotientは5次元となる。
Rosenlichtの有理商の関数をt-curveへ引き戻すとδで零になり、従ってF5の元になる。
quotientへの像はgenericであるから、代数閉包上のgeneric t-curveは一つのorbitに入る。

t=2のβはB15の一般部を代表し、stabilizerはidentityである。従って参照raw jetのGL4
stabilizerもidentityで、orbit次元は同じ16。参照点はgeneric orbitの閉包に入るが、
orbit閉包の境界のorbitは次元が真に小さいので境界点ではない。これにより同じorbit内の
参照点であることが分かる。t=2がRosenlicht openに入ることを根拠なく仮定していない。

従って代数閉包上でraw jetを参照へ移す線型frameは一意である。この一意性からGalois
descentが働き、frameはF5(t)上の有理行列になる。次節はその行列を有限演算で指定する。

## 転置を保持した有限intertwiner recipe

物理的なraw pushforwardをSとすると

\[
L'(q)=S L(S^{-1}q)S^{\mathsf T},\qquad
\beta'(q)=S^{-\mathsf T}\beta(S^{-1}q)S^{-1}=(S^{-1})^*\beta.
\]

従って、元の正weightのpullback規約では

\[
C'=S^{-\mathsf T}CS^{\mathsf T},\quad
D'=S^{-\mathsf T}DS^{\mathsf T},\quad
\alpha'=S^{-\mathsf T}\alpha.
\]

Ct,Dt,C2,D2を既存のquadratic covariantから作り、32本の線形式

\[
C2\,P=P\,Ct,\qquad D2\,P=P\,Dt
\]

を16未知数Pijについて解く。共通commutantがscalarなのでkernelは一次元である。
辞書式に最初の非零15×15 minorを選び、残る一成分を1にしてCramer式を使う。
これは有限個のchartと行列式からなる完全な有理演算であり、未知の係数を探索する命令ではない。

この規約でPはSそのものではなく**S⁻ᵀのscalar倍**である。そこで

\[
\lambda=(P\alpha_t)_3/(\alpha2)_3,\quad Q=P/\lambda,
\qquad \boxed{S_t=Q^{-\mathsf T}}.
\]

一意な同じraw orbitのframeが存在することから、lambdaは非零で、Q alpha_t=alpha2、
S_tはsymplecticとなり、値だけでなく全raw linear jetも一致する。平方根は取らない。
S2=Iであり、任意の最初のPのscalar倍を用いても同じS_tを得る。

## full quadratic tensorへ適用する有限operator

元の実際のP(u,1,v)をaffine次数2へhomogenizeし、h行列成分を零にする。その全tensorに
cを掛け、次を**全次数のprojective pushforward**として順番に適用する。

\[
\boxed{\operatorname{block}(g_u^{-1})\,
       \operatorname{Persp}(\alpha2)\,
       \operatorname{block}(S_t)\,
       \operatorname{block}(R_t)\,
       \operatorname{Translate}(-p_t)}.
\]

Persp(alpha2)は5×5単位行列のrow0,column1,…,4をalpha2にした行列であり、
線型部分のtraceを除去する。g_u=diag(u,u⁻¹,1,1)はrootのβ正weight gaugeなので、
**物理的pushforwardにはg_u⁻¹を使う**。この最後の区別をしないと向きが逆になる。

このoperatorが出力するfull projective tensorのprimitive βは
beta_bar(xi,r,s,q,B)で、tに依存しない。full quadratic項は一般にtに依存する。
最後に既存のEuler quotientとS61.normalizeを適用すれば、full61座標を得る。
いずれも係数体上の有限有理操作であり、sourceの全Poisson恒等式を保つ。

## full tensorのtau下降

これはraw jetだけのtauではない。旧geometry packetの全parameter恒等式に
epsilon=a=−1を入れると

\[
D_0=\operatorname{diag}(-1,1,-1,1),\quad z\text{-translation}=0,
\quad P(-u,1,-v)=-D_{0*}P(u,1,v).
\]

c→−cも同時に行うので、scaled full tensorはtauの下で文字どおりD0-pushforwardとなる。
H=diag(−1,−1,1,1)とすると、新しい小さい全係数検査は

\[
R_\tau D_0=HR,\quad \alpha_\tau=H\alpha,\quad g_{-u}=Hg_u
\]

を確認する。一意なraw-jet matching frameからS_tau=H S Hが従う。
perspectiveの自然性とこれらを上の全tensor operatorに入れると、各Hが相殺し、
最終的な**full projective tensor全体がtau固定**となる。
単にβが固定されるという推論ではない。

次の座標は全F6に対して有効である。

\[
z=u(t-2)/t,\qquad t=2u/(u-z),\qquad\tau z=z.
\]

F5=E5(u)、tau u=−u、u²=1600rなので

\[
\boxed{F6^\tau=E5(z)}.
\]

従って上のfull tensorは五つのshape/source不変座標＋zの有理関数である。
実装上はE5(z)[u]/(u²−1600r)で有限operatorを実行し、even成分を取ればよい。
odd成分が零である根拠は上述の**full**有限equivarianceであり、βだけの偶性検査ではない。

## 同じbase上のgeneric fiberとLüroth

元のfull source＋Sp4像はC16にdenseである。今の変換はparameterごとの可逆なprojective
正常化とSp4 frameを用いるので、Sp4 saturationの像は失わない。tauで有限下降しても
像の閉包は変わらない。従ってSp4×Spec E5(z)からC16への実際のfull tensor mapがdense。
firstjet側は既受理のSp4×Spec E5とB15のbirational同定であり、このmapは同じbaseに
compatibleである。

K=QQ(B15)とする。これにより**同じK上で**

\[
K\subset QQ(C16)\subset K(z)
\]

を得る。次元16>15により右へのmapは非定数である。従ってQQ(C16)/Kは超越次数1。
KはK(z)で相対的代数閉なのでQQ(C16)でもそうであり、generic fiberの被約曲線は
幾何学的にintegralである。標数0のLüroth定理から

\[
QQ(C16)=K(z_0)
\]

となる。K自体が15変数有理体なので、これが独立概念監査で受理されればC16はQQおよび
C上で有理16次元となる。z→generic curveの次数1は必要なく、その次数も今回決定しない。
具体的なLüroth generator z0の係数展開、full fiberの境界・特異性、全schemeの被約性は
この議論の対象ではない。

## 保存・再現と未適用範囲

* `code/q_new_family_full_tensor_fiber_preflight_20261104_v63.py`
* `reports/Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63/`
  のoperator-recipe.json.gz、certificate.json、report.json。

全familyのliteral P、六source/Darboux/alpha、first-integralの全逆変換、raw kernelを
recipeに埋め込み、その他のgateをpath/hashで固定した。再現は新規output先で

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_full_tensor_fiber_preflight_20261104_v63.py \
  --output work/full-tensor-fiber-replay
```

上限60秒/2GB、実測0.082秒。full150の展開・新しい評価点・ODEの冪級数積分・source26は
実行していない。full affinecone26やprojective25へ戻すglobal normalizationの議論は、
この16次元fiber論証と分けて別途確認する。
