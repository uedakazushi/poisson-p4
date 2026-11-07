# L(1,1,2)の実対数源：全unimodular低次延長の非優越

対象は閉じた二次二形式γに対する全延長α=J+β+γで、Jは非退化、
αは閉、Pf αは非零定数である。体積形式で双対なbivectorをP=C+L+Qと書く。
このノートは**実際に対数表示を持つ最高部の全源**を扱う。
対数parameter写像のbase locusで生じる、別表示を必要とする極限は末尾に明記して残す。

## 1. 全二次式の五型への還元

独立な二一次式をx,yとして

\[
 \theta=xyq\left(a\frac{dx}{x}+b\frac{dy}{y}+c\frac{dq}{q}\right),
 \qquad a+b+2c=0
\]

と置く。qは任意の斉次二次式であり、既約性・非退化性は仮定しない。
s=c−a、t=c−bとすると、定数倍を除いて

\[
 \gamma=(s-t)q\,dx\wedge dy+s y\,dx\wedge dq+t x\,dy\wedge dq. \tag{1}
\]

s=0又はt=0なら一次Casimirを持ち、s=tなら `γ=s d(xy)∧dq` は二次pencilである。
これらはそれぞれ `UNIMODULAR_LINEAR_HAMILTONIAN_SUPPORT.md`、
`UNIMODULAR_REGULAR_QUADRATIC_PENCIL_SUPPORT.md` で全延長を処理済み。
以下ではtを1へ規格化してu=s/t、u≠0,1とする。

qの二次元部分空間x=y=0への制限の階数で分ける。
階数二なら、z,wの線形変更とx,yを加えるshearで `q=zw+f₂(x,y)` とできる。
この場合、全bivector最高部は

\[
 Q_{12}=0,\ Q_{13}=b xz,\ Q_{14}=-b xw,\quad
 Q_{23}=c yz,\ Q_{24}=-c yw,\ Q_{34}=(b+c)zw+F_2(x,y),
 \quad bc(b+c)\ne0 .                                   \tag{2}
\]

`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` 第12節の五辺証明は任意のF₂を加えても成立する。
実際、三次Pfaffianで平方を含むz,wの単項式は、L34以外の第三変数係数を全て消す。
残る十二係数へのdiv Lと三次Pfaffianの二式は同じで、L12=0かつ
`L=dQ(a)+(p x+r y)∂z∧∂w` となる。平行移動後はL²=0、
二次PfaffianはC34以外のCを全て消すので、定数Pfaffianも零である。
従って制限階数二の全延長は空。

制限階数一では、平方完成で
`q=z²/2+w m(x,y)+f₂(x,y)` とできる。
m=0なら三変数pullbackで、入口第16節の全支持定理を適用する。
m≠0はx,yの交換・個別scaleを除いてm=x又はm=x+y。
wへx,yの一次式を加えれば、f₂のmで割れる部分を消せる。
商の残る一係数を零又は一へscaleして次の四型を得る。

制限階数零では `q=z m(x,y)+w n(x,y)+f₂(x,y)`。
m,nのspanが高々一なら三変数pullback。
独立ならz,wを変えてm=x,n=yとし、shearでf₂を消す。
従って未処理の全実源は次の五型で尽くされる。

| 型 | q |
|---|---|
| A0 | z²/2+xw |
| A1 | z²/2+xw+y² |
| T0 | z²/2+(x+y)w |
| T1 | z²/2+(x+y)w+x² |
| M | xz+yw |

最初の二一次式が従属する場合はθを直接微分するとγが共通一次因子を持つ。
その全延長は入口第15節で処理済み。q=0も同じ既済零最高部に入る。
これで実対数表示の全配置・全residueを保持した。

## 2. 全30低次係数の有限証明書

公開runner `code/classification_unimodular_logarithmic.py --complete` は、各型で
全24一次係数へ `div L=0` と全三次Pfaffian係数を課す。
平行移動の全imageを計算し、非零minorで選んだ真の線形条件によってのみこれを除く。
平行移動で定数項も変わるが、元Pfaffianは非零定数なので全基点で非退化性を保つ。
その後も全六定数係数を保持し、二次Pfaffianから解けるものだけを消す。
全残余の二次・一次Pfaffianと定数Pfaffianを公開する。

最初の全係数体計算は `CLASSIFICATION_UNIMODULAR_LOGARITHMIC_20261107.json` に保存。
そのproducer bytesを保持したまま、normalized kernelと全例外residueを含む証明書
`CLASSIFICATION_UNIMODULAR_LOGARITHMIC_COMPLETE_20261107.json` を作成した（1.777秒）。
有限契約は `jobs/unimodular_logarithmic_complete.json`。
基礎体QQ(u)の全kernelと全minorを使い、minorや分母の全零点をQQ上で再処理する。
u=0,1とt=0は第1節の既済枝なので、この表から省略してよい。

| 型 | 一般uの一次行列rank | 必要な例外residue | 結論 |
|---|---:|---|---|
| A0 | 18 | 2,5 | 全点で定数Pfaffian恒等零 |
| A1 | 18 | 1/3,2 | 同上 |
| T0 | 18 | −1,2 | 同上 |
| T1 | 18 | 2,3 | 同上 |
| M | 16 | 1/2,3/2,2 | 第3節で全primitive評価rank≤3 |

一般uという記述は任意標本を意味しない。QQ(u)上の全係数恒等式であり、
その一つの非零pivot minorと全表示分母が指定する開集合を覆う。
全例外は表の十一点で、計十六の全線形源を保持している。
A0,u=5では一次行列rank17へ落ちるが、増えた全sourceを残しても定数Pfaffianは零。
A0/A1,u=2では平行移動のimageが三次元へ落ちるため、正常sourceは三変数となる。
これらを通常の四次元translationとして商にしていない。

## 3. M型の全延長は一次評価が非単射

この型の最高部は

\[
 Q=(0,-xy,x^2,u y^2,-u xy,(u-2)xz+(2u-1)yw)             \tag{3}
\]

である（順序12,13,14,23,24,34）。
公開証明書の正常変数をn0,n1,n2,n3、定数自由変数をa0,a1とする。
u∉{0,1,1/2,3/2,2}では、全残式に

\[
 \frac{u-1}{u(2u-3)}n_0^2,\quad
 \frac{3u-2}{u(u-1)}n_1^2,\quad
 -\frac2u n_1n_3
\]

があり、さらにn0=0の下で非零定数倍のn1n2が残る。
従ってn0=0、u≠2/3ならn1=0。
u=2/3でn1≠0を仮定すればn2=n3=0となり、
定数Pfaffianはa0²/u、残る一次係数は
`a0 n1(4u−3)/(u(u−1))`。
非退化性からa0≠0なので矛盾し、ここでもn1=0である。
このu=2/3はQQ(u)表示の全分母が非零の点であり、新しい標本推論ではない。

u=1/2,2の全例外源にも非零定数倍のn0²,n1²があり、同じ結論を得る。
u=3/2ではn1²,n0n2,n0n3が残式となる。
n1=0の後で、一次式はa0n0/4、定数Pfaffianは
`2a0(a0+2n2n3)/3` だから、非退化性よりn0=0も従う。
これで全residueにおいてn0=n1=0を証明した。

この代入後の一般uの全Lは

\[
 (0,0,\tfrac{n_2}{u-1}x,n_3y,
   -\tfrac{u n_2}{u-1}y-n_3x,n_2z+n_3w).              \tag{4}
\]

u=3/2では(4)へ直接代入でき、公開例外sourceと一致する。
u=1/2,2では追加の自由係数が最後の34成分のy又はxに現れるだけである。
従って全例外を含め、z,wの二列は同じ一方向∂3∧∂4へ入る。
元に戻す平行移動微分dQ(a)も(3)から同じ性質を持つ。
よって元の全一次評価L:V→Λ²Vの階数は高々三。
体積形式による同型はこの階数を保ち、unimodular性からLは既にprimitiveである。
入口第14節のinjectiveなW16の例により、これは真の閉条件である。

以上から、**実際にL(1,1,2)対数表示を持つ最高部の全unimodular延長はGを優越しない。**

## 4. 独立検証と残る境界

統合担当は五型への全幾何還元、制限rank二の一様排除、M型の全支持論理を独立に読解確認した。
`CLASSIFICATION_UNIMODULAR_LOGARITHMIC_REVIEW_20261107.json` は
producerをimportしないFraction(u)実装で全十六源を再構成した（2.405秒）。
全行列、全kernelとrank minor、translation全像と真slice、六定数消去、全残Pfaffian、
全分母の十一例外、およびM型のn0=n1=0後の全translationを含む評価rank≤3を独立確認した。

対数型parameter写像には、例えばq=x²,u=3やq=xy,u=1でbase locusがある。
そこで原式が零となった後の再scale極限は、同じ実表示の素朴な代入だけでは得られない。
本証明はこの未知の境界を除去したと主張しない。
外部六成分分類との接続、および既済R13/R22/S23の全閉包は入口第16節を参照する。

## 5. 対数型閉包を覆うquartic integrating factorの全境界

以下は全境界の次の入口であり、この節だけで全低次延長を排除したとはしない。
θの三次係数と、四次多項式Dについて

\[
 D\,d\theta-dD\wedge\theta=0                         \tag{5}
\]

は斉次閉条件である。Dが二一次式と一二次式の積となる射影集合は、
三つのprojective空間の積からの、basepointを持たない積写像の像なので閉じる。
このD集合とθの空間とのincidenceを射影すれば、再び閉集合となる。
L(1,1,1,1)とL(1,1,2)の全generic源が(5)を満たすから、
**両成分の全閉包は、ある非零Dについて(5)を満たす。**
元の対数表示が零となるbase locusでも、この主張は変わらない。

用いる部分分数補題は次である。複素有理閉一形式ωの極が
`D=∏ f_i^{m_i}`だけにあるなら、

\[
 \omega=\sum_i\lambda_i\frac{df_i}{f_i}+dF,
 \qquad\lambda_i\in\mathbf C,\qquad
 \operatorname{pole}_{f_i}F\le m_i-1 .                \tag{6}
\]

これはgenericな一変数について部分分数分解し、閉性から単純極のresidueが定数となることを使って証明できる。
単純極を引いた有理関数をその変数で積分すると有理原始関数を得る。
根を分裂させるために有限代数拡大を用いた場合は、正規化traceで原始関数を元の体へ下降する。
差はその変数に依存しない閉形式なので、変数の数について帰納する。
最終的なdFは元の極以外の極を持たないためFにも余剰の極はない。
f_iに沿うFの極次数k>0ならdFの極次数は正確にk+1であり、掲示上界を得る。
ωがEuler-horizontalかつweight零なら、Euler微分からd(EF)=0。
有理関数の一つのscale軌道上の微分が非零定数になることはないのでEF=0であり、
Fはweight零、`Σ λ_i deg f_i=0`となる。
統合担当はこの部分分数・trace下降・帰納、極次数上界とEuler帰着を独立に読解確認した。

ここでω=θ/Dに(6)を適用する。Dが被約なら、全源は第1–4節と入口第12・16節の
実対数表示に入る。非被約Dの全型と必要な追加源は次の通り。

| D | θ/Dの全源（定数は無視） |
|---|---|
| x⁴ | d(f₃/x³) |
| x³y | λ(dx/x−dy/y)+d(f₂/x²) |
| x²y² | λ(dx/x−dy/y)+d(f₂/(xy)) |
| x²yz | λ dx/x+μ dy/y+ν dz/z+d(f₁/x), λ+μ+ν=0 |
| x²q, q既約二次 | λ dx/x+μ dq/q+d(f₁/x), λ+2μ=0 |

f_kは全k次斉次多項式。x,y,zは表のDで相異なる一次因子であり、線形独立は仮定しない。
最初のx⁴は `γ=4 dx∧df₃` なので既済R13へ入る。
x³yでは(1)のs=3,t=1へ `−4λx² dx∧dy` を加えた全源、
x²y²では `2d(xy)∧df₂−4λxy dx∧dy` を得る。
これらのλ=0は既済実対数源であり、非零λの全源が次の有限判定対象となる。
最後の二型と、これら新源の未処理部分は全て保持する。

## 6. x³yとx²y²の全境界源

以下の全source算術は `CLASSIFICATION_UNIMODULAR_LOGARITHMIC_BOUNDARY_20261107.json` に保存した。
同じrunnerの `--boundary`、契約 `jobs/unimodular_logarithmic_boundary.json` で再現できる（1.174秒）。
`CLASSIFICATION_UNIMODULAR_LOGARITHMIC_BOUNDARY_REVIEW_20261107.json` は
producerをimportしないFraction(u)実装で全十源の行列、kernel、translation真slice、全六定数、
全残Pfaffianと全分母根0,−1/2の再処理を独立確認した（1.644秒）。

非零magnetic係数は全bivectorのscaleとqの変更により一へ規格化できる。
x³yでは(1)のu=3へx² dx∧dyを、x²y²ではu=1へxy dx∧dyを加えたものを使う。

qの制限rank二で、x³yは第1節の任意F₂を許した五辺排除をそのまま適用できる。
x²y²のこの枝は
`Q13=xz,Q14=−xw,Q23=−yz,Q24=yw,Q34=F₂(x,y)`、Q12=0となり、
F₂のxy係数eは非零である。
線形場z∂z−w∂wの重みで全Lを分ける。
weight±2項は三次Pfaffianのxz²,yz²,xw²,yw²係数から全て零となる。
translation以外のweight+1項は

\[
 L_{12}=A z,\quad L_{14}=B x+A F_{yy}y,\quad
 L_{24}=A F_{xx}x+(Ae-B)y,\quad L_{34}=-Ae z,
\]

という形である（Fxx,FyyはF₂のx²,y²係数）。
weight−1も同様でL12=A' w。
二次Pfaffianのz²,w²係数は非零定数倍のeA²,e(A')²なので、A=A'=0となる。
weight0の全解はtranslationと `(p x+r y)∂34`。
translation後の二次条件C∧Q=0はcross四成分を消し、さらにC12 F₂=0からC12=0。
よって定数Pfaffianは零であり、この全枝も空である。

制限rank一・零の還元は第1節と同じで、magnetic項を保つscaleを用いる。
x³yのm=x+yではf₂の残るx²項をu=3の演算子が消すため、T0一型だけが必要。
x²y²のm=x+yでは残るx²係数をmagnetic係数と同時にscaleできないため、
`q=z²/2+(x+y)w+u x²` の全一変数を保持する。
他のaligned型は残るy²係数を零又は一へ規格化できる。
mixed rank零の補集合は全pullbackであり、mixed rank二はq=xz+ywとなる。

| D | q型 | 全残余支持からの結論 |
|---|---|---|
| x³y | A0,A1,T0 | 定数Pfaffian恒等零 |
| x³y | M | n0²=n1²=0、二列比例で評価rank≤3 |
| x²y² | A1 | n0²=0、その代入後2n2³=0、定数Pfaffian零 |
| x²y² | T(u)、u=0,−1/2も含む | (n0+n3)²=0、その代入後2n2³=0、定数Pfaffian零 |
| x²y² | M | n2²=(n0+n3)²=0、非退化性から全正常変数零、評価rank≤3 |
| x²y² | A0 | 下記の実六変数全源と全域次元上界 |

x³yのM型ではn0=n1=0後、L13=0、L14,L23,L24はx,yだけの式、
L34=n2z+n3wとなる。Qのcross成分もx,yだけなので、全translationを戻してもz,wの二列は比例する。
x²y²のM型ではn2=0,n3=−n0後の全条件に

\[
 n_0(a_0+n_0n_1)=n_1(a_0+n_0n_1)=0,\qquad
 \operatorname{Pf}C=a_0(a_0+n_0n_1)\ne0
\]

があり、全正常変数が零となる。同じ二列比例で全一次像は非単射。

最後のA0型では、全残条件と非退化性は正確に

\[
 n_0=n_3=a_0=0,\qquad n_1n_2\ne0,
 \qquad\operatorname{Pf}P=n_1n_2^3
\]

となる。従って固定した最高部

\[
 Q=(0,-x^2,xz,xy,-yz,-xw+xy)                          \tag{7}
\]

の全非退化延長は二正常変数と全四translationの像で、次元は高々六。
ここでは一次評価は基点移動により単射にもなり、rank≤3の主張は使わない。
代わりに、非零線形場 `N=x∂z−z∂w` はx,yとq=z²/2+xwを保ち、div N=0なのでQを保つ。
scalar Euler場もQを保つため、GL4でのQのstabilizerは少なくとも二次元、軌道次元は高々十四。
最高部の非零scaleをさらに一変数として戻しても、全affine飽和支持の次元は
`14+1+6=21` 以下である。
非退化な定数bivectorの六次元GL4軌道への評価は同変で、全点がその軌道へ写る。
従って定数部を任意の固定Darboux Cに固定した支持は次元高々十五となり、
primitive β像も高々十五である。
これは一点での微分rank不足ではなく、全源・全基点を保持した次元上界である。

統合担当はA0のstabilizer、全源・scaleを含む次元上界、全基点での六次元評価差引を独立に読解確認した。
以上によりx³yとx²y²の全境界源もGを優越しない。

## 7. x²yzの全源：全係数の独立再構成完了

第5節のx,y,zとf₁が従属なら三変数pullbackへ入る。
独立ならf₁=wと置き、微分した最高部は

\[
 Q=(X+aN)\wedge(Y+bN),\qquad
 X=x\partial_x-2y\partial_y+w\partial_w,\quad
 Y=y\partial_y-z\partial_z,\quad N=x\partial_w,
\]

となる。ここでa=2μ−λ、b=ν−μであり、λ+μ+ν=0の全平面から(a,b)の全平面へ写る。
三場は可換でtrace零。
wのscaleでNの係数を一斉にscaleできるため、b≠0ならb=1,a任意、
b=0,a≠0ならa=1に規格化できる。a=b=0は既済の五辺型。

全二源の元24一次係数と六定数係数の有限計算は
`CLASSIFICATION_UNIMODULAR_LOGARITHMIC_DOUBLE_LINES_20261107.json` に保存した（0.206秒）。
同じrunnerの `--double-two-lines`、契約 `jobs/unimodular_logarithmic_double_lines.json` による。
b=1の全係数体QQ(a)では一次行列rank20、非零pivot8で、全Lはtranslationのみ。
定数行列はrank6、選択minor−2でC=0となる。いずれのminorにもaの例外はない。
b=0,a=1では一次行列rank18、pivot−4で、translationを除いた全源は

\[
 L=(n_0x+n_1w)\partial_y\wedge\partial_z,\qquad
 C=a_0\partial_y\wedge\partial_z.
\]

両源の定数Pfaffianは恒等零である。
`CLASSIFICATION_UNIMODULAR_LOGARITHMIC_DOUBLE_LINES_REVIEW_20261107.json` は
producerをimportせずFraction(u)と疎多項式で両源を完全再構成した。
全kernel・rank minor・translation真slice・六定数・全Pfaffianが一致し、
全表示分母はaについて非零定数だった。従って例外値なしにこの全型は空である。

## 8. x²qの全支持と対数閉包の非優越

最後のquartic integrating factor型D=x²qも、全係数と全補集合を保持して閉じた。
幾何的還元・全源表・独立算術は `UNIMODULAR_LOGARITHMIC_POLE_BOUNDARY_REVIEW.md` 第3–5節にまとめる。
f₁がxの倍数、μ=0、restriction rank一以下のpullbackを既済の全源定理に接続し、
残る微分二形式を

    γ=2q dx∧dy+(y+x) dx∧dq−x dy∧dq

へ規格化する。Jordan場x∂x+(y+x)∂yのcentralizerを保つ変更で、
restriction rank二はzw+y²+u x²、zw+xy、zw+x²、zwの四型、
rank一はz²/2+xw、z²/2+xw+y²、z²/2+yw、z²/2+yw+x²の四型、
rank零はxz+ywへ尽くされる。任意のbinary quadraticと全mixed係数を残している。

producerと独立Fraction(u) checkerは、全九型と唯一のu=0例外、計十源の
全24一次係数、完全kernel、translation真slice、六定数係数、全Pfaffian残式を再構成した。
最初の九源は全Pfaffianが恒等零。最後のmixed源は全残式からn₀=n₁=0を得る。
その後、元の全translationを戻してもz,wの二つの一次評価列は同じ一方向へ入るのでrank≤3。
従ってこの全源もprimitive β像を優越しない。

以上と第5節の閉incidence・部分分数補題を結合すると、
実L(1,1,1,1)、L(1,1,2)の**全閉包**の全unimodular延長はGを優越しない。
これは未知の全rank4成分を尽くしたとの主張ではない。
E(3)主軌道の全延長空性は `UNIMODULAR_EXCEPTIONAL_SUPPORT.md` で独立検証済みだが、
その実軌道境界からの全接続は別の義務として残る。
