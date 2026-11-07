# rank 二 quadric cone：全 tangent source と有限 Ore 入口

2026-11-07。主義務 G、producer `plane_collision_review`。
実既約被約 quintic D の非正規曲面成分がquadric cone
\(S=(h=q=0)\)、\(q=u_0u_2-u_1^2\)、vertex coordinate vとする。
S上のrestriction rank二を扱い、零divisor Rの被約性・既約性を仮定しない。
本稿は全sourceの幾何入口であり、有限operatorだけでG完了とはしない。

## 1. coneでも全 tangent Euler liftを取れる

Poisson invarianceからraw contractionsはS上でλE、μEとなり、
\(\lambda\in H^0(O_S(1))\)、\(\mu\in H^0(O_S(2))\)。
smooth quadricと同じnormal-line twistsである。
vertexではdq0なのでμ(vertex)=0。
従ってμはv2項を持たないquadraticとしてlift出来る。
qの三つのgradientはu0,u1,u2を張るので、全そのようなμはW(q)としてlinear Wで表せる。
W(h)=λも独立に指定出来、両contractionをEuler gaugeで同時に消せる。

S_reg上ではrank二とactual PfaffianのI2 membershipからconormal bracketは零。
SはCIなのでI/I2はvertexでもlocally free、Sはintegralだから全Sへ零として延びる。
global tangent liftについてstrict \(\{h,q\}\in I^2\) が成立する。
degree二のh-rowを \(h\ell_i+q a_i\) と書くと、modhのconormal条件は
\(\sum a_iq_i=0\)。
従ってaはqのradical、vertex方向だけに残る。
a0はactual hfactorで排除され、a≠0の全sourceだけを残す。
smooth quadricでa全零とした段階をconeへ拡大していない。

## 2. 全 restriction九係数を Jacobian bivector で保持する

cone Sはnormal Gorensteinで \(\omega_S=O_S(-2)\)。
S_reg上のbivectorは \(O_S(2)\) のsectionとしてvertexを越えて延びる。
全restrictionはdimension九のquadratic F(u,v) modqで表せる。
四変数(u0,u1,u2,v)上のJacobian bivector

\[
 J(q,F)^{ij}=\operatorname{sgn}(i,j,k,l)(q_kF_l-q_lF_k)
\]

（k,lは残る二index）を取ると、h-componentは零、qとFはstrict Casimir。
\(E_4\) のcontractionとEuler formulaにより、S_reg上のprojective restrictionは
非零constant倍のFに対応する。
従って全九係数が保持される。
rank二restrictionの条件Fmodq非零を課す前の大きいsourceを用いても上界には安全である。

任意のtangentπからこのJacobian liftを引くとS上零となる。
既受理の全rankzero Koszul/Euler sourceを使い、元の全πは

\[
 \Pi=J(q,F)+hL+qA \tag{1}
\]

で表せる。Lは全quadratic raw sourceのh-coefficient、全50係数、
Aは全constant skew5×5、全10係数、Fmodqは全9係数。
raw source69、Euler冗長五、projective source64である。
restriction mapによりJacobian九方向はrankzero sourceから独立なので、全Euler overlapも既存五方向のまま。
S_regでrestriction零のambient sectionはS全体で零になるため、vertexにだけ支持する余分なsectionを落としていない。

## 3. 全 conormal 条件と actual A 正規化

J(q,F)はh,qの両方にstrict Casimirなので、全conormal linear条件は既rankzero sourceと同じ十一式。
raw69源ではそのsolution dimensionは58、Euler五を引いて53。
Aのh-rowで許されるのは \(A^{h,v}=a\) だけである。
a0はHfactorで排除済み。
a≠0では既受理の実shear
\(u_i'=u_i-b_i h/a\)、inverse \(u_i=u_i'+b_i h/a\)、
さらにvertex scalingでa一とA^{u_i,v}零へ直せる。
qの変化はhlinearとh2だけで全hLへ吸収され、Jacobian restrictionの九係数は可逆に動くだけ。
unknown higher termsを固定していない。

残るconstant partは
\(A=\partial_h\wedge\partial_v+C\)、Cは三非radical変数のconstant bivector。
SO(q)とtensor/residual scalingの全inverseによりCは
\(C=\partial_0\wedge\partial_2\)（nonzero discriminant）、
\(C=\partial_0\wedge\partial_1\)（nonzero nilpotent）、又はC0の全三型。
**rank二ではC0を排除しない**。旧rankzeroのC0 actualfactorはF九係数があると流用出来ない。
この全normalizationは(1)の全F空間を保つ。

## 4. 必要な finite linear→original-Jacobi job

Cを上の三whole型に固定し、L50とF9全59係数を残す。
全350 projective Jacobiをmodhへ取ると、これらのunknownにaffine-linearとなる。
hLとhLのbracketはhで割れ、J(q,F)とhLのcrossもmodhで零、
J(q,F)自身は全FについてJacobian Poissonなので、未知同士のquadratic項が無いためである。
旧conormal十一式も加え、Euler五をexplicit quotientして全affine sourceを保持する。

この全sourceをchart h1へ写し、Ore形
\(P=B(u)+(Y(u)+vT(u))\wedge\partial_v\) が得られるかを全係数で検証する。
手Jac(h,ui,uj)はFのv2係数を零にし、Jac(h,ui,v)はLの該当係数をF_vとCから決める。
Fのuだけの五成分も残り、旧rankzero25変数Ore sourceへ単純にrank証明を流用しない。
全sourceのoriginal350 Jacobi・全126actual-Pfaffianを復元し、
homogeneous BのPoisson条件で得る追加の全linear constraintsを確認する。
これが次の一回有限jobの範囲で、33変数等のGB又はslice反復は行わない。

## 5. 将来の全G上界へ戻すnormalizer次元

Aut(S)の実dimensionは十二、tensor scalarを含め十三。
allowed A-stratumのaffine dimensionsはC semisimple七、nilpotent六、zero四。
従ってnormalized A stabilizerのdimensionsは六、七、九。
各normalized full Poisson source dimensionをdとすれば、全PGL image上界はそれぞれ
\(18+d,17+d,15+d\)。
rankzero packetと同じactual A-normalizationのfiber countであり、F九係数の可逆移動を含む。
したがってd≤6,7,9が最弱の十分条件。
finite affine sourceやnormal symbol因子だけからこれを満たすとは主張せず、
full Jacobi/whole Ore supportを別に証明する。
Rのreduced/nonreduced、C0、F_v零、F0零等の全境界を保持する。
