# nil31短translationの全H+既約quartic源はMPS：独立検証

2026-11-07。課題Eの型4、\(Z=y\partial_x+z\partial_y+\partial_w\)。
新しいfocused reviewer codeだけを作り、producer・CAS・legacy算術をimportせず、
既存の独立Fraction/Laurent多項式と五変数疎算術で再構成した。

現行codeは classification_nil31_short_quartic_review.py、
現行証明書は CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_REVIEW_V2_20261107.json。
有限契約は同prefixの REVIEW_JOB_20261107.json である。
既存/shared/global/Gitには手を加えていない。
最終実行約1.074秒、最大RSS20,832 KiB、CPU上限170秒、address space上限1,536 MiB。

## 1. 公開inputとimmutableな初回証明書

入力は公開済みの NIL31_SHORT_QUARTIC_SOURCES, SUPPORT, AMBIENT, COMPONENT,
RECOGNITION の五report、および元のNILPOTENT_SOURCESとNIL31_SHORT_SUPPORTである。
初回reviewはprivateな完成packetを読んだためV1をimmutableに保存した。
公開RECOGNITIONとprivate完成packetについて、全structured recognition結果と
producer source SHA256

    7f4559458fb49cb34e49d0032d3e6b35a582444df1bdeed417f8d80289bd6ce4

が完全に同じであることを確認し、V2では公開reportのみへread/bindingを移した。
private work fileは現行code/jobの依存ではない。
既済の全semisimple型12支持定理のnote/reportと、型12のスペクトルが
\((-2,-1,1,2)\) である全線形source表もinput hashへ追加した。

## 2. 全九coreから十quartic候補への完全還元

元の110個のdivergence/Lie線形行から全90係数を再び解き、
六括弧・四Jacobi・実Pfaffianを再構成した。
全Jacobi係数のmonic集合とcore reportの全13式が一致する。
既済の独立九core theoremと、そのGroebner不要の
\(\tau=ct,\ J=cj,\ K=c^2k\) 恒等式を再利用する。

\(c=t=k=j=0\) の二つの連続coreでは、全Jacobiから

\[
 \begin{array}{ll}
 b=1/3:&(3a+1)n_8/6=0,\\
 a=-b:&(1-3b)n_8/9=0
 \end{array}
\]

を得る。\(n_8=0\) のとき最初のmagnetic式は
\((a+1/4)n_7=0\)、最高四次Pfaffianは \(a n_7z^4\)。
従って四次性が残り得るのは

\[
 (a,b)=(-1/3,1/3),\quad(-1/4,1/3),\quad(-1/4,1/4)
\]

だけである。これらの恒等式を元の全source上で検査した。
非零cの五枝、非零tの一枝、非零kの一枝と合わせ、全十候補となる。

各候補でlower六変数とmagnetic六変数の全二段行列・右辺を元Jacobiから別に作った。
選択したpivot minorは各caseの非零core parameterの単項式と非零有理数だけであり、
指定localizationの全点でunitである。
nonpivot変数が自由のまま残ること、全pivot解が全元Jacobiへ戻ること、
残るcompatibilityが保存された集合と正確に一致することを検査した。
特に候補7の \(n_{13}n_8=0\) は捨てていない。

## 3. 九候補の普遍除外と唯一の完全六源

候補0,1,2,4,5,7,8,9では、実Pfaffianの全parameter恒等式として、それぞれ

    z, z², 3cz+1, 15cz−1, z+3t, z², z, z⁴

の非定数因子を持つ。全quotientを展開して検査した。
係数の特殊値でquarticの次数が落ちる場合も、単純H＋既約quarticの対象には入らない。
候補6は最高四次Pfaffianが恒等零である。

唯一残る候補3は

\[
 c=n_{17}\ne0,\quad a=-5/12,\quad b=-1/12,\quad
 t=-1/(72c),\quad j=7/(72c),\quad k=-1/(216c^2)
\]

であり、自由変数は
\((c,n_{12},n_{13},n_6,n_7,n_8)\) の正確に六個。
全十八置換、六括弧、実Pfaffian、全Jacobi恒等零を確認した。
parameter空間は全 \(\mathbb G_m\times\mathbb A^5\) である。

実Pfaffianはxに一次で、その係数は

\[
 c\left(z-\frac1{6c}\right)^3.
\]

この係数の唯一の根 \(z=1/(6c)\) での残りは

\[
 -\frac{(12n_{12}c+2n_{13}+6cy-w)^2}{3456c^3}.
\]

w二次係数は常に非零なので、xの二係数のgcdは一である。
Gaussの補題により、**全parameterで**実quarticは絶対既約。
\(xz^3\) の係数cが非零なので全体次数も常に四であり、斉次化でh因子は増えない。
一般parameterだけの既約性ではない。

## 4. 正規点から出るliteralな全Poisson曲線

\(c=1\)、他の五自由変数が零の点から、元の全六式を432倍してPを作り直した。
保存された真normal R、補正Tの全150係数とも一致を確認した。
新しい形式的modelを仮定せず

\[
 P(e)=P+eR+e^2T
\]

そのものについて、全350個の射影Jacobiと全十個のraw斉次Jacobiが
多項式として恒等零になることを別実装で確認した。
全eでh-Casimirも保つ。

全homogeneous cofactorを元のdivergenceから再構成すると

\[
 (0,\ 432y,\ 432z-1296ex,\ e(216y+36w),\ 432h-5184ex).
\]

その全5×5行列から特性多項式を別に展開して

\[
 u(u^2+93312e)(u^2+373248e)
\]

を得た。\(e\ne0\) ではaffine線形部が可逆でtranslationを消せ、
スペクトル比は \((-2,-1,1,2)\)、既済のregular semisimple型12となる。

実C-star pushforward恒等式

\[
 \operatorname{diag}(1,s^2,s,1,s)_*P(e)
       =sP(e/s^2)
\]

も全monomialで検査した。
tensorのpushforward weightはP,R,Tでそれぞれ \(1,-1,-3\) であり、
逆のweight conventionに依存せずこのliteral identityを検証している。
非零eは全てPGL/scalar同値となるが、以下の閉性認識には一般eだけで十分である。

## 5. MPSへ入る最弱の橋：全源は正規点のorbit closure

P(0)のquarticが既約でh因子を持たないので、一般eでもその開条件が保たれる。
従って一般P(e)は、既済の型12の**全支持定理**により有限閉MPS unionへ入る。
これはMPSの一意性や、SNC族の一覧だけからの外挿ではない。
有限閉unionの閉性からP(0)もそこへ入る。

元の全20×20 affine cofactor作用を再構成すると、Zのorbit次元は13。
完全六源とhyperplane選択を動かした射影飽和閉包の次元は高々

\[
 6+13+4-1=22.
\]

この全閉包は、既約な \(\mathbb G_m\times\mathbb A^5\) と連結PGL5からの像の閉包なので既約。
一方、P(0)の実PGL/gauge/scale行列はrank47であり、
vertical rank24とscale一を除く実射影orbit次元は22である。
従ってこの22次元orbitの閉包は、全六源のPGL飽和閉包と一致する。

MPS unionはPGL安定な閉集合でP(0)を含むので、
そのorbit closure、従って全六源とそのPGL飽和を含む。
これで型4の全H+既約quartic支持のMPS被覆が得られる。
全源を一つの実orbitだと主張する必要はなく、orbit closureで十分である。

## 6. 追加で確認した全周囲のsmooth23証明

上の最弱橋には不要だが、producerの強いrank証明も独立に確認した。
全150方向から全350×150 Jacobi微分Dを生成し、
非零102minorをFractionで再計算した。
全GL5二十五方向、vertical二十五方向、scalar一方向を作り、rank47、
Rを加えた実際の多項式曲線方向を含む行列のrank48を確認した。
全350×48のD積は零なのでDのrankは高々102でもあり、正確に102である。

projective tangent次元は \(150-24-102-1=23\)。
PGL/scalarとliteral curveが与える実familyの接rankから次元下界も23となる。
従ってPは周囲Poisson空間の滑らかな23次元点で、通る全成分は唯一である。
この唯一性と既約六源を用いる別の帰属橋も成立するが、
第5節のorbit closureによる弱い証明を主として採用する。

higher jet、nilpotence exponent、formal IFT modelの全形式完全性は使用しない。
新成分を追加する結論でもない。
