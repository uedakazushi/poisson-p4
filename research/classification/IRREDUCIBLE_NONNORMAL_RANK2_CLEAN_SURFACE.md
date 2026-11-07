# rank二の滑らかな二重曲面：vertical成分なしの全枝をprojected Veroneseへ還元

2026-11-07、主義務G。既約・被約・非正規五次退化因子Dの非正規曲面成分Sについて、
次の全subbranchを扱う。

1. Sは滑らかで、その一般点でπのrankは二。
2. DのS一般点での横断型はordinary node。
3. `b:Y=Bl_S P4→P4`、exceptional `p:E=P(N_S/P4)→S` とすると、
   strict transform D' とEの交叉divisorはSを支配する成分だけを持つ。
   曲線へ落ちるvertical divisorは無い。

正規化X又は導手Cが滑らかであるとは仮定しない。
Sが特異、非node横断型、vertical交叉成分がある場合は補集合として残す。
結論は、この全subbranchのSが一つのsmooth projected Veronese軌道に入るという有限還元である。
その全Poisson源は独立有限計算でaffine次元≤2、全PGL射影像≤22と証明された。
従ってこの全subbranchは未知W優越成分の一般支持にはならない。

## 1. blowupへのPoisson liftとexceptional residue

SはHamiltonian derivationで不変なので、ideal I_SはPoisson idealである。
Sの一般点で接方向Poissonが非退化で、実PfaffianはI_S²に入る。
局所normal座標u,vと接座標s,tを取ると、mixed成分はI_Sに入り、
Pfの一次normal項は接方向の非零係数とnormal-normal成分の一次項の積である。
従って一般点で `{I_S,I_S}⊂I_S²`。
これはconormal bundle上のO_S双線形bracketであり、滑らかな既約S上で一般点で零なので全域で零。

chart v=uτで `{u,τ}={u,v}/u`、
`{s,τ}=({s,v}−τ{s,u})/u` はいずれも正則である。
従ってπはYへ正則にliftし、Jacobiもdense openから全域へ延びる。
`−K_Y=5H−E`、DのSに沿う重複度は二なので、liftの実退化因子はD'+Eである。
その逆形式のE上residueをβとすると、βはE上の閉有理一形式で、
divisorialな極は `Γ=D'|E` に限られる。

ordinary nodeの二normal枝のbiresidueをλとする。
一般点のlog skew行列による既存局所計算から、π|Sのrank二はλ≠0と同値。
blowup後の二交点でβのresidueはλと−λである。
閉性により各既約polar divisor上のresidueは定数である。
従って二交点が同じ既約bisectionに属することは出来ず、Γのhorizontal部分は
相異なる二つのdegree一成分Γ₁,Γ₂へsplitする。
仮定3とgeneric nodeからΓ=Γ₁+Γ₂で、これ以外のdivisorial poleは無い。

## 2. Gysin関係と次数四の強制

滑らかなprojective E上の閉有理一形式のresidue関係は

`λ[Γ₁]−λ[Γ₂]=0 in H²(E,C)`。

Γᵢが特異でも、polar divisorのresolution上で通常のlog residue関係を取り、
exceptional divisorを押し下げれば同じ関係を得る。
Γᵢ自身のsmoothness又はSNCを追加仮定する必要は無い。

`ξ=c₁(O_E(1))` とする。`E|E=−ξ` なので

`[Γ₁]+[Γ₂]=2ξ+5p*H`。

Γᵢのfiber degreeは一、Picard群のprojective bundle公式より
`O_E(Γᵢ)=O_E(1)⊗p*Bᵢ` と整数line bundle Bᵢを取れる。
従って `B₁+B₂=5H`、`c₁(B₁)=c₁(B₂)`。
`A=B₁−2H` と置けば

`2A≡H`、Aはample、`4A²=H²=deg S`。

既存の全normalization-length公式からdeg S≤6。
A²は正整数なので、deg S=4、A²=1となる。
ここではcohomology classの一致を無断でlinear equivalenceに読み替えていない。
数値同値と整数intersectionだけで十分である。

## 3. smooth degree四とH≡2AからP2を初等的に得る

SがP3に入ればsmooth quartic hypersurfaceでK_S=0。
A²=1のRiemann--Rochは `χ(A)−χ(O_S)=1/2` となって矛盾する。
より小さなspanはdegree四のsurfaceではない。
従ってSはP4で非退化である。

一般hyperplane切断curveのdegreeは四、P3で非退化なのでCastelnuovoのboundから
その種数g_Hは0又は1。随伴公式は `g_H=3+A K_S`。
AのRiemann--RochからA K_Sは奇数なので、g_Hは偶数、従ってg_H=0、A K_S=−3。
hyperplane restrictionのH¹単射からq(S)=0。
K_S A<0なのでp_g(S)=0、従ってχ(O_S)=1。
Riemann--Rochとh²(A)=h⁰(K_S−A)=0からh⁰(A)≥3。

|A|に非零fixed divisorがあれば、ample Aとのintersectionで
A²≥1+1となるので不可能である。
basepointがあればA²=1から一つのsimple basepointだけで、blowup後のmobile systemの平方は零。
これはcurveへ写り、そのimageは少なくともP2で非退化なのでdegree≥2。
一般memberが少なくとも二つの有効curveの和となり、再びample AとのintersectionでA²≥2となる。
従って|A|はbasepoint freeである。

Aはampleなのでその射は有限、A²=1なのでimageのdegreeと射のdegreeは共に一。
imageはP2、有限双有理射は同型である。従って `(S,A)=(P2,O(1))`。
Pic(P2)=Zより、ここで初めて `H=O_P2(2)` とlinear equivalenceも得る。
P4への五切断はVeronese P2⊂P5の一点射影であり、smooth embedding条件から中心はrank三。
従って既存と同じ唯一のsmooth projected Veronese軌道に入る。
Fujitaその他の未確認のpolarized threefold分類は使わない。

## 4. rank零の証明を転用せず、完全なtangent sourceへ

rank零PVは既済だが、ここではπ|Sは一般rank二なのでその排除は適用しない。
新しい全条件は、既存parametrizationq(a,b,c)に対して

`q∧Π(q)=r₃(a,b,c) (dq_a∧dq_b∧dq_c)`。

r₃の全十係数を未知数に加え、ambient全150係数と合わせて全280×160行列を作った。
rank121、kernel39で、Πへの射影は単射、Euler24をquotientした完全tangent sourceは15次元。
全350 Jacobiのparameter-quadratic row rankは54。
実anticanonical quinticの全係数はこのquadratic rowspaceには入らないので、
この判定からrank四排除は主張しない。

証明書は `CLASSIFICATION_PROJECTED_VERONESE_TANGENT_20261107.json`、
runnerは `code/classification_projected_veronese_tangent.py`、実行前gateは
`jobs/projected_veronese_tangent.json`。全sourceのproducer bytesを保存した。
固定PV全Poisson源のaffine次元は、全54 Jacobi quadricが生成する理想からの
242polynomialの全10,079 QQ lift恒等式と、leading monomialのminimum cover13により≤2である。
完全GBの再実行又はrank零への制限を使わず、元Jacobi理想への各membershipから直接上界が従う。
PVを保つPGL群はprojection centerのconicを保つPGL3、すなわちidentity成分SO3の三次元群。
PGL×projective sourceへの自由fiber作用により、固定源次元dなら全像次元は高々20+dとなる。
従って全像次元は≤22となり、この全subbranchの未知優越成分を排除した。

dimension証明書は `CLASSIFICATION_PROJECTED_VERONESE_TANGENT_DIMENSION_20261107.json`。
独立 `code/classification_projected_veronese_tangent_review.py` は、
全280×160、kernel39、Euler24商15、全350 Jacobi、実quintic、全10,079 liftと
全242leading monomialをFractionだけで再構成してPASS（11.508秒、最大RSS59,680KiB）。
独立証明書は `CLASSIFICATION_PROJECTED_VERONESE_TANGENT_REVIEW_20261107.json`。
dimensionと独立reviewの実行前gateは同prefixの`jobs/projected_veronese_tangent_{dimension,review}.json`。

## 5. 残す補集合

本稿は非正規曲面全体のsmoothness、ordinary node性又はvertical divisorの不存在を証明していない。
特異S、非node横断型、Γのvertical成分、rank零曲面は別に残る。
正規化X・conductorの全点smoothnessも主張しない。
第1–3節はwhole subbranchから唯一PV軌道への証明であり、
第4節の固定源dimension/supportの独立算術と合わせ、明記した全subbranchのG用閉鎖が成立する。
