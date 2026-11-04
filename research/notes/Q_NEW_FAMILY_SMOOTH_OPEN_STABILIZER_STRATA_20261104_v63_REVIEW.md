# uvp開集合全域のstabilizer：独立監査

2026-11-04。原著`Q_NEW_FAMILY_SMOOTH_OPEN_STABILIZER_STRATA_20261104_v63`を独立に監査する。
領域はu v(32u+v)≠0。原著producerをimport・再実行せず、新しい点・素数・GBを使わない。
独立checkerはFraction係数の独自8変数Laurent多項式演算とAST parserだけを使い、
原著のnested fraction-field polynomial ringと別に全homogeneous係数を作る。

## intrinsicな旗とconstant compatible空間の跳び

保存済みの元のfull affine PからSchouten bracketを直接微分し、constant bivectorの
20×6行列の全120成分を作り直す。原著行列との完全一致を確認する。
v,2,−64,4/3を係数に持つ四行だけで

\[
k_{xy}=k_{xw}=k_{yw}=0,\quad k_{xz}=-u k_{zw}
\]

が必要となる。そのvectorを全20行へ戻すと、残る唯一の条件は
(4/3)(u²−64)kzw=0である。kyz列は零。
従ってconstant compatible空間はu=±8で二次元、それ以外で一次元。
跳びで増えるbasisは−u∂x∧∂z+∂z∧∂wである。この例外を無視したCの一意性は使用しない。

元のPのdivergenceからhomogeneous Nを独立に再構築する。像はey,ez平面に入り、
rows(y,z),columns(h,x)のminorは4u²/3なので、全u≠0でImN=span(ey,ez)。
full stabilizerのaffine行列Mはdivergence covarianceからNM=κMNを満たすため、この像を保つ。
唯一のtriple HとQ|Hの唯一のvertexも全uvp領域でintrinsicなので、元の一般affine旗は
compatible Cが跳ぶu=±8を含めて正当である。

κを含む一般旗からκf=1、b=64c、d=a、D=a²−64c²、κ=1/D等へ落とす既受理の
係数比較はu,v,Dの可逆性だけを使っている。genericなu²−256やv+16uの反転をここへ
紛れ込ませない。

## 全quadratic caseの独立確認

Q1,Q2,Q3の全多項式を独立エンジンで処理する。通常のcombination
Q1+Q3+(v−p)Q2=−32ac(u²−256)を確認する。
u=16s（s=±1）についても、全係数で

\[
Q_2=(a+8sc)^2(1-(a-8sc)^2),
\]
\[
Q_3+vQ_2=16c(256+sv)(a+8sc)
\]

を確認する。D≠0ならa±8scは共に単元なので、v+16u≠0ではc=0,a=±1しか残らない。

J=0、すなわちv=−16uでu²≠256ならac=0であり、a=0の場合は
Q2=64c²(1−64c²)だからc=±1/8。この二familyが新しい候補となる。

共鳴点u=16s,v=−256sではa−8sc=epsilon=±1、a+8sc=t∈Gmなので
a=(t+epsilon)/2、c=s(t−epsilon)/16となり、これが全quadratic候補である。
Q1,Q3のQ2倍数性も独立に照合し、他のcaseを落としていない。

## 全150係数、unit系、共鳴Gmの棄却

新しいJ=0二familyと共鳴四caseについて、元のfull Pをhomogenizeし、独自Laurent環で

\[
\kappa M\,P(X)M^{\mathsf T}-P(MX)
\]

の全10bivector成分×15homogeneous monomialを直接作る。
J=0二familyでは全150係数がQQ(u,T)上で零となる。有限点検査ではない。

残るH,K,Zのconstant系は原著のpulse法を再実行せず、P(0)とtarget translationの値から
直接導く。初めの二行のH,K部分は

\[
\begin{pmatrix}-512\kappa c/9&8\kappa a/9\\
8\kappa a/9&-8\kappa c/9\end{pmatrix},
\]

第三行のZ係数は4u/3。他の二つの第三行成分もconstant bivector P(0)から直接作る。
全matrixとrhsが保存値に一致し、detは−256u/(243D)、全表示領域で単元である。
保存されたH,K,Zを代入した全三線形式も一致する。
従って特殊parameterで未確認のrank低下を仮定していない。

四つの共鳴Gm候補で実際に残るfull150残差を作り直し、原著の全六非零係数と一致させる。
各caseで全残差はt²−1の倍数であり、少なくとも一係数は±(t²−1)/24というunit倍。
従って**全tensorを保つための追加条件は正確にt²=1**であり、二次条件だけからのGmは残らない。
これは別の係数点を試す操作ではなく、QQ(t,T)上の全多項式の比較である。

## groupとexhaustion

保存Sを独立に取り直し、S²=E²=I、SE=ES、SN=−NS、EN=NEを確認する。
さらにJ=0の二つの全T-familyについて

\[
M_+(T)=\exp\bigl(3TN/(4u)\bigr)S,
\]

と、M−がESにmodular Gaを掛けたfamilyであることを、matrix全成分で照合する。
M−のtimeは3(T−(ES)yh)/(4u)であり、平行移動のoriginの違いも省略しない。

旧a=±1 familyは全QQ(u,v,T)恒等式として受理済みで、分母にu²−256やv+16uを持たない。
今回の全域unit constant系が一意性を保証するので、その二familyは全uvp域で同じく尽くす。
J=0では新しい二familyを加え、共鳴点でもt=±1から同じ四blockへ戻る。

以上から、全u v(32u+v)≠0について原著の次の分類を**独立に受理する**。

| 領域 | actual algebraic projective stabilizer |
|---|---|
| v+16u≠0 | Ga × μ2 |
| v=−16u、u²=256も含む | (Ga ⋊ μ2) × μ2（最初のμ2がGaを反転） |

componentのlinear blockは互いに異なり、Gaと余分な有限factorを同一視しない。
標数0のstabilizer groupはsmoothなので、これはgeometric pointsだけの数え上げではない。
QQ上の各stratumでもsplitし、Cおよび幾何的base extensionで同じ結論となる。
u=0、v=0、p=0の境界、parameter間のKlein four同定とは区別する。

再現は同名独立checkerに新規output先を渡す。原著のsource/evidenceは変更しない。
独立証拠は全新caseの150残差、unit系、flag修復、有限group恒等式を保存する。
