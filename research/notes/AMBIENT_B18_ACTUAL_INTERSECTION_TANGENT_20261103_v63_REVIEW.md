# B18原点のactual相対接空間：独立150×39監査

2026-11-03。`AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63` のCA/CB/D接空間と
dim/e/depth上限を承認する。v62によるorbit depth下限は別監査であり、ここでは繰り返さない。

## 元のRでのcotangent分解

受理済みの旧共通target R29上で `A=R directsum I*t`、t²=nu、I⊂m_R、nu∈m_R。
`m_A²=m_R² directsum m_R*I*t` なのでedim A=29+mu_R(I)。
intrinsicなedim39は新しい全形式モデルで既に確定しており、この等式は旧Rにも適用できる。
従ってmu_oldR(I)=10。新モデルの座標を旧CA/CB相対baseへそのまま移していない。

既受理のseven linear・two WU/WV・one U³はJ/mJで独立であり、mI⊂mJなのでI/mIでも独立。
mu=10によって全basisとなる。従って全actual odd coefficientのbasis展開の定数は、
全firstjet L、zero-firstjet quadratic Q、monic critical U-axis cubic C3で正確に決まる。
critical U曲線では全odd generatorがO(U³)なので、mIの像はO(U⁴)となりcubic評価に混入しない。
以前に未証明だったI∩mJ=mIを先に仮定する必要はない。

## even embeddingとsource kernelの二乗

生CA/CBの20basisをcanonical化・h/x交換して再構成する。
P0の全25 GL列は、外積indexと二次monomialの逆座標微分から独立に作る。
旧受理済みの19orbit＋9fiber列J28を全150係数で照合する。

純粋kernel curveのfiber係数と、有限group `I+pM+p²M²/2` を作り直す。
CBではratio×fの二次項も保持する。
独自の疎係数環とliteral inverse-coordinate/exterior pushforwardによって、
全150の一次係数零と二次Wを確認する。producerのSymPy有限pushは再利用しない。
Rはactual target関数なのでR(p)=O(p²)、trace-zero t(p)=O(p)、I⊂m_RからtI=O(p³)。
従ってこのWはgenuineなeven一次embeddingの一方向である。
一般のsource二次jetを新しいtarget tangentとして採用したものではない。

共通29row blockのdetはCA9216/25、CB−576/25。
J28とWからeven embedding E150×29を独立に復元する。
odd embedding B150×10と合わせ、selected row blockの三角構造を利用して39×39 inverseを使わず
正規化する。新producerのE/B/Tと、共通正規化150×39行列の**全係数**が一致する。
selected39rowがidentityなのでrank39、CAとCBの同一行列性から和も交わりも同じ39平面。
actual ideal sumに対して `T_D=T_CA intersect T_CB` なのでT_Dも39である。

## 次元・重複度・depthと範囲

既知のreduced supportはincomparableなN28とR20（cone）。従ってDのdim28。
topdimはNだけで、そのgeneric actual length1とP0でのreduced Nのmultiplicity2からe(D)=2。
低次元Rやembedded成分の長をtop multiplicityの和へ加えない。
minimal prime R20もAssに入るのでdepth-localization inequalityよりdepth D≤20。
scalar smooth factorを除くprojectiveはdim27/edim38/e2/depth≤19。

この接空間結果からD_redの接空間、nilradicalのlinear visibility、full原点環、exact depthは推論しない。
特にT_N=37だけを見てreduced union N∪Rの接空間も37だとはしていない。
対象は指定P0とそのopen B18 orbitのG移動であり、他の境界点へ拡張しない。

独立コードとhash証拠：`AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63_REVIEW`。
