# 全1+1+3支持から未知W優越成分を排除する

2026-11-07。課題G。実退化因子 `D=h x F₃`、二超平面は相異なりF₃は既約三次式とする。
`TWO_HYPERPLANES_CUBIC_SUPPORT.md` の全normal模型と全nilpotent配置を入力に、
既存MPSへの帰属又は次元24以下の有限全源によって、未知W優越成分の一般支持を排除する。
Eの残る全源帰属は別に残し、この結果から新成分を数えない。

## 1. 全支持の有限な残余

旧全模型は `P=x∂x∧V+W`, `V=By+cx+d`、Wの全45二次係数を保持する。
Bの半単純部が非零の全支持は、旧ノート第2–4節によりMPS40の閉和に入る。
Bがnilpotentの場合、同第5.1節の全十一正常形を使う。

Jordan3の三型とJordan21の `(L,M)=(0,0),(h,0),(h+x,0)` はF₃可約又はJacobi不整合で空。
`(0,h)` のu依存枝もF₃がhで割れるため空。
u非依存の枝は全平面十三配置に入る。
`(h+x,h)` のp=±1全二源は既知の同じ22次元全成分へ入る。
残る非平面源は次の通りである。

| 残余 | 固定Vでの完全源の次元 | sourceを保つ群の次元 | 全PGL像次元上界 |
|---|---:|---:|---:|
| `(L,M)=(h+x,h)`, p=0 | 5 | 6 | 23 |
| `(L,M)=(h,x)` | 6 | 6 | 24 |
| `(L,M)=(0,h+x)`, n₆=1, 旧(30) | 6 | 7 | 23 |
| B=0 rank二translation、非平面L=h型 | 8 | 10 | 22 |
| 同、非平面L=h−x型 | 7 | 9 | 22 |

最後の二行を落とさないことが重要である。
B=0 rank二translationには実際に非平面支持がある。
例えば元21係数でn₁₁=n₁₇=1,n₄=−2、他零なら
`F₃=h²u−hxv+x³` は既約で、全Jacobiを満たす。
この二行の完全導出は第4–6節に示す。
B=0 rank一translationはF₃が一次式で割れ、rank零ではPのrankが四にならない。

## 2. 最初の三全源の大域的なpreserving group

`V=v∂u+L∂v+M∂w`、L,Mが独立な場合、次の六線形場を取る。

`v∂u+L∂v`, `w∂u+M∂v`, `h∂w`, `x∂w`, `h∂u`, `x∂u`。

全てh,xを点ごとに固定し、Vと可換である。
六場は独立でLie bracketで閉じ、順序(u,v,w,h,x)で全行列はstrict upper triangular。
従って標数零の有限exponentialにより六次元の連結代数的unipotent群Kを作る。
この群は固定Vの全Jacobi支持、二超平面、F₃の既約性を保つ。

`(h+x,h)` の全Jacobi支持はp=0,±1という三つの相異なる閉じた既約五parameter源である。
連結群は有限個の既約成分を置換しないため、p=0全源が保たれる。
`(h,x)` の全既約三次支持は唯一の完全六parameter源に入るため、同じKがその全対象を保つ。
点における群作用rankを全源へ外挿した議論ではない。

旧(30)では `V=v∂u+(h+x)∂w` であり、次の七場を取る。

`h∂u`, `x∂u`, `v∂u`, `w∂u+(h+x)∂v`, `h∂w`, `x∂w`, `v∂w`。

全てVと可換でh,xを固定し、Lie閉、順序(u,w,v,h,x)でstrict upper triangular。
七次元連結unipotent群Kとなる。
determinantは一、uの係数は一で他の座標にu項を作らないため、
実三次式のu係数 `−n₆(h+x)²` を保つ。
従って正規化n₆=1と、旧(30)の**全**六parameter源を保つ。
n₆を正規化する前の非零値も、旧証明の可逆(u,v)同時scaleでこの源へ入る。
n₆=0の枝は別の平面支持として保持する。

## 3. 平面十三配置の全像は次元23以下

`TWO_LINES_PLANE_CUBIC_SUPPORT.md` の全配置・全六magnetic源を入力にする。
三変数baseの選択はGr(3,5)の六次元、各配置のPGL₃軌道は高々八次元、合計十四。
十三配置は有限で、追加の連続modulusを持たない。
固定したDの完全Saito源はanchorのGL₂ frame四次元と、全六magnetic係数で尽くされる。
projective tensor scaleを除いてfiber次元は高々九。
従って全平面支持の全像は次元高々 `14+9=23`。
fiberへのbase一次式の加算は既に全magnetic二次式へ吸収され、別の自由度として足さない。
この次元論は全十三源に適用し、個別の残るE成分帰属を仮定しない。

## 4. B=0 rank二translationの全二十一係数を保つ模型

`V=x∂u+h∂v`, `T=hu−xv` とする。
元の全21係数sourceは正確に、A,B,Cを(h,x,w)の任意の二次式、a,b,cを定数として

`P_xu=x²`, `P_xv=hx`, `P_xw=0`,
`P_uv=A+hu−aT`, `P_uw=B−bT`, `P_vw=C−cT`。                (1)

元parameterとの対応はa=n₄,b=n₁₁,c=n₁₈で、A,B,Cは各六係数を全部残す。
`L=bh−cx` とすると、実三次式は

`F₃=LT+xC−hB`。                                         (2)

h=1で `G=xC−B`, `T=u−xv` とする。
補助座標(x,T,w,v)ではPは、三変数上の `−F₃∂T∧∂w` と
`Q∧∂v`、`Q=(x,A+(1−a)T,−C+cT)` の和となる。
これは方程式を記述する補助座標であり、非線形座標変更をPGL同値として使わない。
元全Jacobiは正確に次の二polynomial恒等式へ同値である。

`I=xL'+cG_w+C_w L=0`,
`II=xG_x+A L−C G_w−(1−a−C_w)G=0`。                     (3)

実際、元の最初の三Jacobiは恒等零、最後は `−(I T+II)`。
新独立Fraction reviewerは元の全21係数でこの全恒等式を再構成した。

L=0なら(2)は(h,x,w)だけの平面型。
L≠0はwの非零scaleとh/xの相対scaleでL=h又はh−xへ正規化できる。
L=xは二超平面の交換でL=hへ入る。
この交換では同時にu,vも交換し、新無限遠に対してEuler外積を引き、全tensorを−1倍すると
再びcanonical `V=x∂u+h∂v` となる。Tは符号を変え、(b,c)は交換される。
従って無根拠にh-Casimir gaugeを維持したと仮定していない。
より具体的にはh'=1上で、旧A,B,Cへ(h,x,w)=(x',1,w')を代入したものを
Ã,B̃,C̃と書くと、新しい係数は
`A'=Ã`, `B'=−C̃−x'w'`, `C'=−B̃−w'`, `a'=−a`, `(b',c')=(c,b)`。
全六括弧に対する実際のlinear交換とEuler補正から得る表示である。

## 5. L=hの全非平面源は八parameter

b=1,c=0とする。(3)のIからC_w=0であり、Cはxの二次式となる。
F₃のh=0制限はCのx²係数をc₂としてc₂x³。
F₃の既約性からc₂≠0（零ならhがF₃を割る）。
Bのxw,w²係数をそれぞれb₄,b₅とする。
IIはAを一意に

`A=C G_w+(1−a)G−xG_x`

と定める。その三次項は正確に

`−c₂(b₄+a+2)x³−2c₂b₅x²w`。

従ってa=−b₄−2,b₅=0で、全sourceは

`C=n₀+n₁x+n₂x²`, n₂≠0、
`B=n₃+n₄x+n₅w+n₆x²+n₇xw`, `a=−n₇−2`

という全八parameterに尽くされる。
全式の除算は非零c₂と非零定数だけで、c₂=0補集合は既約性に矛盾すると処理した。
上記Aを戻せば全Jacobiは全parameterで恒等零。
F₃はuについてprimitive一次式で、斉次係数h²と残りのn₂x³は互いに素なので絶対既約である。

## 6. L=h−xの全源は七parameter

b=c=1とする。Iから

`B=C−xw+D(x)`、Dはxの二次式。

IIは

`A={C G_w+(1−a−C_w)G−xG_x}/(1−x)`                    (4)

と定める。分子のx=1への制限を取ると、w²,w,1の係数から順に

`C_{w²}=0`, `a+C_w|x=w=0+C_{xw}=0`, `D_{x²}=D_const`

が得られる。各係数pivotは非零定数であり、別のparameter反転は無い。
この三条件は分子が1−xで割れるために必要十分である。
分子の全次数は三以下なので、商Aの次数は二以下となり、これ以上のJacobi条件は残らない。
全七parameterは

`a=n₄`,
`C=n₀+n₁x+n₂x²+n₃w+(−n₄−n₃)xw`,
`D=n₅+n₆x+n₅x²`。

F₃はTについて一次であり、x=1では `G=w−D(1)` が零多項式でないため、
係数1−xとGは互いに素である。従って全parameterで絶対既約となる。
式(4)はparameter又は1−xの主開だけの有理表示ではなく、全parameterで正確な多項式商である。

## 7. B=0二全源のpreserving normalizer

次の八独立線形場はV,h,xを点ごとに固定し、Lie閉なunipotent代数を作る。

`h∂u,x∂u,w∂u,h∂v,x∂v,w∂v,h∂w,x∂w`。

全行列は順序(u,v,w,h,x)でstrict upper triangular。
従って八次元連結代数的群Uとなる。
det=1でTの係数も一だから、b,cの正規化を保つ。
さらに次の対角変換を加える。

`D₁(ρ)=diag(1,ρ,ρ,1,1)`、`D₂(λ)=diag(1,1,λ,λ,λ)`。

D₁はVを固定し、(b,c)を(b,c/ρ)へ送る。
D₂はVをλVへ送るが、全tensorをλ⁻¹倍する同じprojective作用によりVと(b,c)を保つ。
従ってL=hの全八sourceは `U⋊(Gm)²`（十次元）、
L=h−xの全七sourceは `U⋊Gm`（九次元）で保たれる。
全sourceの完全性を使うため、個別の接ベクトルや一点のorbitrankからの推論ではない。
独立Fraction計算は全行列rank、Lie閉性、Vの直線保存とnil部分のstrict triangularityを確認した。

## 8. 自由fiber作用からGへ

上のsource SはcanonicalなPx係数を持つので、P(S)内のtensor scaleは既に固定されている。
従って `PGL₅×S→P(全Poisson空間)` のdomain次元は24+dim S。
各preserving群Kはfiberに `(g,P)→(gk⁻¹,kP)` として自由に作用する。
第一factorへの作用が自由なので、Pのstabilizerを余分に足す必要はない。
像の次元上界は `24+dim S−dim K`、第1表の値となる。
非正規化係数を戻す操作は全て既存PGL又はprojective tensor scaleに含まれ、二度数えない。

既知閉成分へ入る枝と、全像次元24以下の有限個の構成可能sourceの閉包とで、
**全1+1+3支持**を覆った。
未知W優越成分の必要次元は25–27なので、その一般元がこの支持に入ることは不可能である。
この結論は全1+1+3のE成分認識を完了したという意味ではない。

## 9. 証明書と独立確認

* `CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NORMALIZER_DIMENSION_REVIEW_20261107.json`：
  最初の三全sourceの群6/6/7次元を独立Fraction行列で再構成。
* `CLASSIFICATION_TWO_HYPERPLANES_CUBIC_RANK2_TRANSLATION_SUPPORT_20261107.json`：
  B=0の全8/7源、全四Jacobi、実Pf/cubic。producer実行0.176秒。
* `CLASSIFICATION_TWO_HYPERPLANES_CUBIC_RANK2_TRANSLATION_REVIEW_20261107.json`：
  元21係数でのI/II、全8/7源のJac/Pf、群10/9の行列を独立Fractionで再構成、PASS0.046秒。

runnerは `code/classification_two_quadrics.py --two-lines-rank2-translation`、
独立runnerは同`_review.py` の `--reduced-divisor-normalizers` と `--two-lines-rank2-translation`。
三つの実行前契約は `jobs/two_hyperplanes_cubic_{normalizer_dimension_review,rank2_translation_support,rank2_translation_review}.json`。
旧reportは変更せず、追加targetの証明書だけを保存した。
source幾何の完全性、平面全支持、自由作用による次元論は本稿と旧ノートの別証明入力である。
