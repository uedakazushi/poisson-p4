# 全non-normal integral CI22のrank二枝：PGL像二十四上界

2026-11-07、主義務G、producer plane_collision_review。
**命題（producer証明、別作者review前）。** integral非退化non-normal CI(2,2) S⊂P4を
非正規曲面supportに持つ、実退化因子Dが既約被約quinticの全rank四Poisson構造の
PGL5飽和像はdimension二十四以下である。
rank二restrictionが対象だが、計算sourceはrank零も安全に含める。
restrictionの零divisor Rの既約性・被約性・vertical型を一切仮定しない。

## 1. 全五surface orbit

全symmetric-pencil canonical theoremからの網羅性・各modelのintegralityは
`NONNORMAL_CI22_RANK_TWO_SOURCE_ENTRANCE.md` §1で証明した。
normalityを仮定してantiK sectionを五係数へ延ばす議論はここでは使わない。

| model | q1 | q2 |
| --- | --- | --- |
| regular22_1 | 2x0x1+2x2x3+x4² | x1²+x3²+x4² |
| nil32 | 2x0x2+x1²+2x3x4 | 2x1x2+x4² |
| singularM1 | x0x1+x3² | x0x2+x4² |
| cone_node | 2x0x1+x2²+x3² | x1²+x2²+2x3² |
| cone_cusp | 2x0x2+x1²+x3² | 2x1x2+x3² |

regular member側は同じ固有値の二blockが共に長さ二以上の場合だけsingular lineを持つ。
size五から(22)+(1)又は(3,2)だけ。
all-singular側のM1+regular2は唯一モデル、common-radical coneのbase C4はpa一のため
node/cuspだけ（複数singular pointとJordan4のline componentを除く）。
全二/三固有値をpencil PGL2で固定できるため、表の全五モデルは一つずつのactual surface orbitである。

## 2. full150から全sourceを作る

全raw homogeneous quadratic coefficientを保持し、各qiについて
`E∧(Π#dqi)=0 modI` を全componentへ課す。
degree四のI=(q1,q2)は全q_i×quadraticのrank二十九、有限QQ rowspaceでexactな剰余を取る。
次にintrinsic local generators q_i/x0²のbracket numerator

`K=x0{q1,q2}−2q2{q1,x0}−2q1{x0,q2}`

をdegree五のI²（全q1²,q1q2,q2²×linear、rank十五）へ入れる。
nonplaneの既約quinticに対するconormal-abelian補題はgeneric Sへ適用出来、
prime CIのI²-primary性によりdense x0 chartから全homogeneous条件へ戻る。
したがってどのgeneric Poisson点もこの全sourceへ入り、特異curve/conductor上の余分なsectionも落としていない。

Euler24が全kernelに入ることを直接確認し、explicit coordinate minorでquotientする。
完全kernel/free-coordinate逆写像と全元150係数、350Jacobi、実126Pfを保存した。
初回linear/leading判定は次の通りであり、全G成功とは扱わず不変に残した。

| model | 全tangent projective源 | conormal後の全源 | 初回coarse上界 |
| --- | ---: | ---: | ---: |
| regular22_1 | 26 | 10 | 26 |
| nil32 | 26 | 12 | 28 |
| singularM1 | 26 | 10 | 26 |
| cone_node | 26 | 14 | 実D0 |
| cone_cusp | 28 | 16 | 32 |

初回source reportの`all_G24_closed=false`は正しい失敗診断として不変保存する。
それをより細かいleading idealで反復せず、元quadricsの全被約支持へ経路を変更した。

## 3. 全被約支持の有限分岐

完全な方程式分岐・全反転と零補集合・自由parameter数は
`NONNORMAL_CI22_REDUCED_SUPPORT_REDUCTION.md` §§1–4を参照。
以下はその帰結と、最後のactualD確認である。

**regular22_1.** (p2−p8)²0からp2=p8=t。
t非零では全supportは二parameter以下。
t0,p1非零では三以下。
t=p1=0では二つの全五parameter linear Poisson族A/Bへ分かれる。
その全元350は恒等零だが、全126 actualPfも双方**恒等零**であった。
従ってrank四・既約D枝は前二つだけに入り、固定S上でdimension三以下。

**nil32.** 全元quadricsから、p10非零、p10零かつp4非零、双方零の三補集合へ分ける。
前二つの必要sourceは各五parameter以下。
双方零では(p3,p6,p8)が三lineに入り、非零coreのsourceは四/三/三以下。
core零では残るf(g+f)=0の二枝とその全交わりを保持し、最大五以下。
全rational mapの逆caseを手で証明し、分母は各caseで非零としたparameterだけである。
これらは必要supersourceであり、全parameter値がPoissonであると要求していない。
残るcompatibilityはsourceを小さくするだけである。

**singularM1.** p8非零では三以下。
p8零ではp5非零をx=p3/p5のgeneric/x2/x3へ分ける。
p5零もp3非零/零、p2非零/零、p7非零/零を尽くし、全supportは最大五以下。

**cone_node.** 全actualD coefficientは元350JacobiのQQ-linear combinationになる。
全三十三の非零D係数のoriginal functionalsを保存し、全Poisson点でactualD0。

**cone_cusp.** 元式p13*p15=0。
両零枝では全actualD0。
片方非零の全二枝はそれぞれ五parameter linear Poisson族へ入り、各全350をliteralに確認した。
全rank四sourceのdimensionは五以下。

## 4. normalizerを正しく戻す

実surface stabilizerの下界は順に三、四、四、五、五である。
regular22_1の三field（centralizer二とpencil lift一）、nil32の既受理四field、
singularM1のrank二Veronese-center stabilizer四、coneのvertex affine変換五を
入口note §3に明記した。
これらはfixed Sの全sourceを保ち、PGL5 frame側の作用は自由である。

個々のparameter wallからそのままgroup dimensionを引くのではなく、
**全有限coverからstabilizer不変な全rank四sourceのdimension上界**dを得た後で、
`24+(d−1)−dimAutS` を使う。
各枝の単独のAut不変性を仮定していない。

| model | 全rank四・既約D源のd上界 | AutS下界 | PGL像上界 |
| --- | ---: | ---: | ---: |
| regular22_1 | 3 | 3 | 23 |
| nil32 | 5 | 4 | 24 |
| singularM1 | 5 | 4 | 24 |
| cone_node | 空 | 5 | 空 |
| cone_cusp | 5 | 5 | 23 |

これで全non-normal CI22のG非優越が従う。
独立受理済みnormal CI22の全rank二空性と合わせれば、全integral非退化CI22 rank二枝は
未知のdimension二十五〜二十七G優越componentの一般元を持たない。
全点のE帰属又は新成分数を決定したとは主張しない。

## 5. 再現可能な証明書

reusable runnerは `code/classification_nonnormal_ci22_tangent.py`。
初回fullsource modeと`--support-completion` modeを持ち、旧main function ASTは不変。

- `CLASSIFICATION_NONNORMAL_CI22_TANGENT_20261107.json`：全線形source/all150/350/126、失敗coarse上界。
- `CLASSIFICATION_NONNORMAL_CI22_SUPPORT_COMPLETION_20261107.json`：全必要rational mapsと全元式残差、
  regular全五族D0、node全元functionals、cusp二wholeopen/sourceと両零D0。

初回1.538秒、support補完2.269秒。後者は全mapの代入だけからcoverageを推論せず、
別の全case逆写像の手証明と組み合わせる。新規性は未監査。
