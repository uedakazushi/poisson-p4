# Universal Sp4曲線のWeyl下降：幾何的種数3、低六不変量体への次数8

2026-11-04。ordered regular-D chartで証明済みの**全generic曲線**を、実際の
Sp4 orbit-invariant fieldへ下降する。結論はk=QQとCの両方で

\[
F=k(W_{16})^{\operatorname{Sp}_4},\qquad
K=k(F_4,p_1,p_2,p_3,p_4)
\]

に対し、F/Kの幾何的generic曲線の種数は3、かつ

\[
\boxed{[F:K(N_{12})]=8}
\]

である。有限な固有値base変更を除くために、同じ幾何的曲線をWeylで新たに商する
Riemann–Hurwitz計算を行うわけではない。NEW hypersurfaceのfinite multisectionの
次数10/20/40、ambient不変体の純QQ有理性、regular invariant環の生成、全域有限性は
この結論に含めない。

## 1. 正確なordered係数体

S=F4/20、P=p4、Delta=x²−y²とする。regular semisimple Dの固有値を
x,−x,y,−yと順序付きで取ると

\[
x^2+y^2=S,\quad x^2y^2=P,\qquad
a=-\frac{p_3+2x^2S}{2x\Delta},\quad k_1=p_1,\ k_2=p_2.
\]

従って

\[
K'=K(x,y)=k(x,y,a,k_1,k_2).
\]

この五変数表示は新しい無関係なbaseへの置換ではない。既存の全pencil計算と
universal rank6は、五つのspectral base変数の独立性と、この有限表示を保証する。
x,yのsigned permutation group W(C2)は位数8であり、その固定体は
k(S,P,p1,p2,p3)=Kである。よってK'/KはgenericにGalois、正確な次数8。
各Weyl操作後のaは同じ固定p3から上式で再計算する。x,yの平方と、後述のb,e gaugeの
平方根を混ぜない。

## 2. 実際のW16からpair・torus sliceへの同じ函数体

Zを実際のquadratic covariant map beta→(C,D)の被約像閉包とする。既存の全係数恒等式、
kernel-lineの非零minor、quadratic componentの非零性から、generic fiberは正確に
{beta,−beta}である。全てQQ上の式とminorであり、characteristic0でgeneric degreeは2。
Sp4の中心−Iがbeta→−betaを実現するので

\[
k(Z)=k(W16)^{\{\pm1\}},\qquad k(Z)^{Sp4}=F.
\]

ここで既に中心の二重性を数え終えた。後のgaugeへもう一つの二次拡大を加えない。

D regularの開で、DをCartanのdiag(x,−x,y,−y)へ移す。二つのdiagonalizationの違いは
normalizer N(T)であり、通常のrational slice argumentにより

\[
k(Z)^{Sp4}=k(Z_{\rm Cartan})^{N(T)}.
\]

この主張はC上だけに依存しない。固有値を付けた後は各eigenspaceのkernel vectorを
そのfield上で取り、±固有値pair間の非零symplectic pairingで一方をrescaleできる。
余分な平方根を取らずにsymplectic eigenframeとなる。標準のsigned permutation
normalizer代表もQQ上の0,±1行列で存在する。

### b=e=1 gaugeをrationalな選択とは呼ばない

T=diag(t,t⁻¹,u,u⁻¹)はb→t²b、e→u²eと作用する。b,eが単元なら有限étale coverで
b=e=1とできるが、その平方根が元のfieldに属すると仮定しない。
gauge sliceをS_gと書くと、T×S_g→Cartan sliceはT[2]=(mu2)²を用いた標準の
finite étale descentであり、函数体に対して

\[
k(Z_{\rm Cartan})^T=k(S_g)^{T[2]}
\]

となる。これはsplit torus上の商のfield descentで、QQでも成立する。
未商のgauge sliceが複数成分なら、ここでk(Sg)は各成分の函数体の積であるtotal rational
algebraと読み、T[2]の不変部分を取る。未商sliceの幾何的既約性を仮定する必要はない。
後述のrelative-sign quotientの正確な表示と全ordered曲線の既約性により、その不変部分は
実際に一つのfield Lとなる。
対角のcentral mu2はpair(C,D)にtrivialであり、effectiveな残りのmu2だけが
(g,h,i,j)→(−g,−h,−i,−j)と作用する。

既受理のLaurent chartはこの相対符号を

\[
R=g/h,\quad T_0=gh,\qquad
g^2=RT_0,\quad h=T_0/g,\quad i=vg/T_0,\quad j=u/g
\]

で**正確にquotient済み**である。右側のu,vはcross product gj,hiで、NEW族の二parameter
u,vとは別物。従ってordered curveのfield Lはまさにk(Z_Cartan)^Tであり、ここへ
相対符号の二次拡大を再追加してはならない。

## 3. 全ordered曲線とsemilinear Weyl作用

独立受理済みのGENUS3_BRANCHは、以下を同時に証明した。

* 全generic (2,2)交点代数はK'(m)上dimension8。
* reciprocal quotientの四次曲線の正規化は幾何的既約・種数1。
* その上のmonic二次coverは四つの無限遠点でのみ分岐し、幾何的既約・種数3。
* 従って全八次代数はfieldであり、実際のm-dominant選択成分がこの**全field**を占める。

この最後の点により、Weylが選択成分を別の未知成分へ送る懸念は解消する。
L/K'はregularな超越次数1の幾何的曲線函数体で、滑らかなproper模型の種数は3、
[L:K'(m)]=8である。

normalizerの実際のpair actionから、W=N(T)/TはT-quotient field Lに作用する。
この作用はK'を保つが、その上では位数8のsigned permutation Galois作用であり、
**semilinear**である。gauge chartやLaurent openが各Weyl代表でregularに保たれる必要はない。
function field上のbirational actionで十分で、全八つのR,T0への展開式も不要である。
この作用はK'上でfaithfulなのでL上でもfaithful。前節のslice同定から

\[
\boxed{F=L^W}.
\]

三つの異なる二重操作を区別する。

1. beta±はSp4中心によってpairへの移行時に消えた。
2. effective relative gauge signはR,T0への移行時に消えた。
3. reciprocal involutionはR,T0の上の別のinvolutionで、その種数1 quotientからの
   **分岐二重coverをLに保持する**。これをさらにgaugeと見なして種数1へ落としてはならない。

非自明なWeyl要素は全てK'に非自明に作用する。K'を固定するreciprocal involutionは
このWeyl作用に含まれない。

## 4. Galois下降は同じ幾何的曲線を保つ

Artinの固定体定理により[L:F]=|W|=8である。L/FのGalois群Wで、K'の全元を
固定する部分群はidentityだけである。Galois correspondenceから直ちに

\[
\boxed{FK'=L}.
\]

さらに[FK':F]=8=[K':K]なのでFとK'はK上linearly disjointであり、F⊗K K'はこのLとなる。
これは不変式の点評価からfield degreeを推測する議論ではない。

F内の元alphaがK上代数的なら、K'上も代数的である。L/K'はregularだからalpha∈K'。
しかしalphaはW固定なのでalpha∈K'^W=Kとなる。従ってKはFで相対的代数閉。
標数0なのでF/Kはregularで、実際のgeneric spectral curveは幾何的にintegralである。

F/Kの滑らかなproper曲線へK'をbase changeするとL/K'の曲線になる。
幾何的種数はconstant base changeで変わらないから、

\[
\boxed{g_{\rm geom}(F/K)=3}.
\]

Wは一つの固定された幾何的fiberに作用する群ではなく、固有値baseのGalois群として
係数を動かしている。従ってこれはWで曲線を新たに商して種数を下げる操作ではなく、
同じ幾何的曲線のdescentである。

## 5. actual low-six fieldへの次数8

既存のintrinsic formulaは

\[
m=\frac{N_{12}}{2\delta_8},\qquad
\delta_8=(F_4/20)^2-4p_4=p_0\in K^\times
\]

を与える。従ってmはFに属し、K(m)=K(N12)。universal rank6によりmはK上超越的である。
K'(m)/K(m)は同じconstant extensionで次数8。受理済みordered次数と合わせると

\[
[L:K(m)]=[L:K'(m)]\,[K'(m):K(m)]=8\cdot8=64.
\]

一方[L:F]=8なので

\[
\boxed{[F:K(N_{12})]=[F:K(m)]=64/8=8}.
\]

ここで低六不変量はF4,p1,p2,p3,p4,N12である。そのsubfieldはactual invariant field全体
ではなく、正確に指数8の部分体である。対応するrational mapがgeneric degree8であることを
示すが、全域でfiniteなmorphism、全special fiberの八点、regular ring生成を主張しない。
NEW multisectionの既知10/20/40とは対象が異なり、この8をNEWへ読み替えない。

## 保存と検証の範囲

* `code/q_firstjet_universal_weyl_descent_20261104_v63.py`
* `reports/Q_FIRSTJET_UNIVERSAL_WEYL_DESCENT_20261104_v63/`

producerは既存の独立六gateと全dependency/evidence hashを確認し、八つのsigned permutationを
有限に列挙して証明schemaを保存するだけである。CAS、Weyl curve式の新展開、点、素数、
新しい不変式jobは実行しない。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_firstjet_universal_weyl_descent_20261104_v63.py \
  --output work/universal-weyl-descent-replay
```

一次保存は専用の独立概念監査待ちとする。QQ/Cのgenus・degreeを証明しても、
ambient quotientが純QQ有理であるという別の主張を追加しない。
