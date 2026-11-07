# 非退化degree三曲面上の全rank零sourceと、rank四に十分な単一局所化

2026-11-07、主義務G。ここではintegral・非退化なS⊂P4のdegreeが三の場合だけを扱う。
本noteの全sourceとorbit reductionは証明済みである。
以下に明記する二次rowspace testは負であり、それだけからrank四不存在又はG排除を主張しない。

## 1. 二つの曲面軌道と全三十係数

外部入力はminimal degree classificationである。
D. Eisenbud–M. Green–K. Hulek–S. Popescu,
*Small schemes and varieties of minimal degree*, Amer. J. Math.128(6)(2006),1363–1389,
arXiv:math/0404517v5, 著者公開PDF
<https://eisenbud.github.io/papers/pdfs/2006-002.pdf> のTheorem0.1、preprint p.2を
2026-11-07に直接読んだ。冒頭はvarietyをreduced irreducible schemeと定義し、
smoothnessを仮定せず、minimal degreeの全varietyを分類する。
同箇所の古典分類の参照先はD. Eisenbud–J. Harris,
*On varieties of minimal degree (a centennial account)*,
Proc. Sympos. Pure Math.46, Part1(1987),3–13である。

degree三・dimension二・ambient P4にこの分類を適用すると、
Veronese及びそのconeはdegree/ambient dimensionが異なり、quadric/linear spaceも該当しない。
rational normal surface scrollの整数型は `a+b=3, 0≤a≤b` なので
`S(1,2)` と `S(0,3)` の二つだけである。各型は一つのprojective orbitである。
その全quadric idealは次の通りである。

| 型 | q0,q1,q2 |
|---|---|
| smooth cubic scroll S(1,2) | x0x3−x1x2, x0x4−x1x3, x2x4−x3² |
| twisted cubic cone S(0,3)、vertex e4 | x0x2−x1², x0x3−x1x2, x1x3−x2² |

これらはそれぞれ行列
`[[x0,x2,x3],[x1,x3,x4]]` と `[[x0,x1,x2],[x1,x2,x3]]`
の三つの2×2 minorである。
height二のHilbert–Burch列
`0→O(-3)²→O(-2)³→I_S→0`
がexactで、syzygy columnsは順に
`(x3,-x2,x0),(x4,-x3,x1)` と
`(x2,-x1,x0),(x3,-x2,x1)` である。
line bundle cohomologyから `h¹(O_S)=0`、`h⁰(I_S(2))=3`、
`h¹(I_S(1))=0` を得る。

従ってEuler補題 `RANK_ZERO_REGULAR_LINEARNORMAL_SURFACES.md` により、
S上でprojective bivectorが零となる**全section**は

`Π=q0 A0+q1 A1+q2 A2`、各Agは任意のconstant skew5×5 matrix

で尽くされる。parameter数は全三十であり、Euler商に残るkernelは無い。

## 2. 全Jacobi二次rowspace testの正確な負結果

`code/classification_cubic_scroll_zero.py` は全三十係数を保持し、
strict JacobiをEulerとwedgeした全350 projective式を、
全465 parameter二次monomialの整数行列として構成する。
実anticanonical式は `D=Σ_i(-1)^i x_i Pf(Π_without_i)` から直接構成し、
全126 coordinate五次monomialの係数を調べた。
`reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_JACOBI_20261107.json` に
全行列、全exact RREF、元の350行による各RREF functional、
各D係数のfunctionalと非零residualを保存した。

| 型 | 全Jacobi row rank | Dの全係数のrowspace包含 |
|---|---:|---|
| S(1,2) | 99 | false |
| S(0,3) | 88 | false |

これは線形membershipという一つの十分条件が成立しないという結果である。
非零residualはPoisson点ではなく、rank四点又は未知成分の存在を意味しない。
rank二tensor `q*A`（A decomposable）の族は全固定source内に実際に残るため、
全rankを含むsource次元だけでは必要なG排除に足りない可能性がある。

## 3. smooth scroll：dimension六の群と一つのrank四patch

scrollのparametrizationを
`[u s,u t,v s²,v st,v t²]` とする。
base PGL2のdimension三、二fiber blockの相対scalingのdimension一、
`u→u+v(αs+βt)` の二unipotent parameterは、scrollを保存する忠実な群Hを与える。
projective次元は `dim H=6` である。
最後の変換はambient座標で

`(x0,x1)→(x0+αx2+βx3, x1+αx3+βx4)`

となり、`(x2,x3,x4)` を変えない。
`q2=x2x4−x3²≠0` のとき、二unipotent parameterでx0=x1=0へ動かせる。
その後base PGL2は非退化binary quadraticのprojective空間に推移的に作用する。
従って `{q2≠0}` はP4内の単一のdense open H-orbitである。
代表を `p0=e3` とすると、q0=q1=0、q2=−1であり、

`D(p0)=−Pf((A2)_without_3)`。

全rank四tensorのDは非零五次式なので、このdense open orbit全体に零とはならない。
従って、どのrank四tensorもSを保つHの元によって `D(p0)≠0` のpatchへ動く。
初めの座標でD(p0)=0の補集合を捨ててはいない。

## 4. cone：dimension八の群と一つのrank四patch

最初の四座標はtwisted cubicのSym³C²表現で、base PGL2が作用する。
さらにvertex方向の相対scaling一つと
`x4→a x4+b0x0+b1x1+b2x2+b3x3` の四translationがある。
これらはconeを保存するdimension八の群Hを与える
（GL2のdimension四、vertex scaling一、translation四からprojective scalar一を引く）。

最初の四座標が相異なる三根を持つbinary cubicを表す点集合はdense openである。
PGL2の三重点推移性とvertex translationにより、これは単一H-orbitとなる。
代表 `p0=(1,0,0,1,0)` ではbinary cubicはs³+t³、q0=q2=0、q1=1なので、

`D(p0)=Pf((A1)_without_0)−Pf((A1)_without_3)`。

従って§3と同様、cone上の全rank四sourceもこの単一 `D(p0)≠0` patchのH飽和で尽くされる。

## 5. 局所化dimensionがGへ与える最弱の結論

上の各patchを、全350 Jacobi式と `inverse*D(p0)−1=0` によるaffine scheme Xとする。
inverseはD(p0)から一意に定まるため、d=dim Xはrank四open source自体の次元である。
`PGL5×X→P H⁰(E)` の各fiberには、
`(h,λ)→(g h⁻¹,λ hπ)` によるdimension `dim H+1` の部分がある。
patchを保つhの集合はHの非空openなので、全Hがpatchを集合として固定する必要は無い。
従ってprojective像の次元は `24+d−1−dim H` 以下である。

smooth scrollでd≤7、coneでd≤9を証明すれば、各全rank四飽和像は24以下となり、
未知25〜27次元W優越成分の一般元から排除される。
これは全rank四不存在、全E帰属、又は各固定源の完全成分分解より弱い十分な目標である。
局所化に対する一回ずつの有限native計算をjob
`jobs/cubic_scroll_zero_rank4_dimension.json` で定義し、
exact liftとleading monomialのdimension上界だけをsuccess criterionとした。
有限制限に達した場合には、このdimension結論を主張せず、全sourceとpatchを再開状態として保存する。

## 6. 有限計算の結果と、停止した代数路線

`CLASSIFICATION_CUBIC_SCROLL_ZERO_RANK4_DIMENSION_20261107.json` に、
各orbitのnative入力と結果を保存した。
どちらも一回のnative liftstdが60秒制限に達し、dimension証明を得なかった。
この失敗をG非優越の証明とはしていない。

局所化による全次数の混在を避けるため、別の一回だけの有限目標として
`D(p0)²∈J` をtorusのdegree四weight blockで検討した。
これは指数を求める研究ではなく、上の全rank四open-orbitを閉じる一つのscalar radical十分条件である。
smooth scrollでは唯一のweight blockが152×240、rank142であり、平方targetはそのrowspaceに入らない。
全exact residualは `CLASSIFICATION_CUBIC_SCROLL_ZERO_SQUARE_SMOOTH_SCROLL_12_20261107.json` に保存した。
coneは三weight blockのうち最初の1092列のblockで100秒制限に達した。
`CLASSIFICATION_CUBIC_SCROLL_ZERO_SQUARE_TWISTED_CUBIC_CONE_03_RUNTIME_20261107.json` は
producer bytes、input hash、command、SIGALRM終了142を保存し、数学的結論が無いと明記する。

ここでこの代数路線を止め、高い冪、長いjet又は同じ全ringのGBを反復しない。
両degree三rank零曲面のG排除は未完了であり、§3–5の完全なpatch還元を保持して
normal-cone又はvertical residueの幾何へ移る必要がある。
