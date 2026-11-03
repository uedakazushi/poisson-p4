# Q_b第一CB枝とCAの相対CI4・実quotient diagramの独立監査

2026-11-03。対象はclosed transverse CI4とrelative CI4の二packet。
98秒のSymPy workerは再実行しない。既受理の独立QQ(b)係数-array evaluatorへ
今回の保存二normal frameを代入し、元の全150係数から直接二次jetを再計算する。
元の二次差のsource方向をそのまま使う計算ではない。

CBの22 elementary basisをraw mu*f/Hからcanonical化し直し、CA20も生のCAから再構成する。
両sourceの全29方向、58frame、36pivot inverse、20baseと二normalのkernel条件を確認する。
normal二方向とその和で、実5×5群と逆座標代入のt²係数を求め、polarizationによって
全150×3 raw quadraticを得る。CBのratio×binary mixed項はformulaに保持する。
今回の両normal ratio speedが零であることは別に確認する。
36implicit補正後の全150×3を保存generic値に照合し、全150行が二つの選択row63/68の
QQ(b)線形結合になることを完全係数等式で検証する。

四つの三次multipleのdetは `1024 b^4/[25(b²−3)^4]`。
全frame/inverse/raw/corrected arrayのb=1特殊化もpoint保存値へ結びつける。
pointはQ1=6(r−z)²、Q2=−23r²/15+18rz/5−9z²/5、三次span det64/25。
従ってclosed二式のidealは `(r²,(r−z)²)`、HF(1,2,1)、length4。

## actual idealと相対自由性

20baseは実smooth Sの座標、二normal idealをmとする。
全S上の接次元22とsectionのstraighteningから、全residualはm²に入る。
QQ(b)上のfull quadratic rank2はGL-equivarianceとdominant G×b像によりS上の全3次minorを零にする。
b=1の二次minorがunitなので、完成base R20上でも全quadratic rowは選択二式のR-linear span。
これは一つのmodpあるいはpoint上のrankから上限を外挿する議論ではない。

三次spanのunitから、actual Fi=Qi+O(m³)に対しm³⊂(Fi)+m⁴。
有限S-module `m³/(m³ intersect (Fi))` へNakayamaを使うとm³⊂(Fi)。
従ってQiも(Fi)へ入り、逆にQiの三次multipleも全m³を生成するので、actual ideal=(Q1,Q2)。
全残り150residualの二次項をこの二式で除けばm³となるため、他のhigher residualも全て入る。
graded piecesはR,R²,R,0,...であり、全相対ringはR-free rank4。
非退化なmultiplication pairingのunitとcomplete char0 base上のHensel/GL2によって
`R[u,v]/(u²,v²)` へ正規化できる。

## quotientの実際のbinding

元の二CB meeting方向を新しいframeのnormal二方向へ投影する。
独立コードはphysical tangentから20actual base方向を引いた全150残差をQQ(b)で確認し、
normal像が `c0=−b³+b,c1=−2b²` であることを結びつける。
従ってkernel lineは `eta=c1*r−c0*z` であり、任意の抽象length2 quotientを選んだものではない。
eta²と `(b*r−z)²` の選択二quadricによる表示を厳密に確認する。
変換detと後者のmeeting tangent値は `−b²(b²−3)`、b=1でunit。

quotient先はdouble smooth Sなので、Rのparameterをその像で取ればR[ρ]/ρ²である。
normal変数の像は必ずc_i(R)ρで、二次以上は零。
kernel lineのisotropic性はintrinsicなquadratic diagramの性質であり、QQ(b)恒等式と
GL-equivariance/dominanceで相対全体へ延びる。
unitな二isotropic lineを選び、quotientをu->0,v->ρとして

    R[u,v]/(u²,v²) -> R[v]/v²

を得る。kernel(u)=Ru directsum RuvはR-free rank2でsquare-zero。

この結果はCB第一枝とCAのintersectionおよび実meetingへのquotient diagram。
全actual Dへの非flat base change、Tor、nilmodule、primary分解、旧P_aのlabel交換は一切承認しない。
b²=3は選択frameの分母であり、幾何的例外を分類したという意味でもない。
保存先は `reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63_REVIEW/`。
