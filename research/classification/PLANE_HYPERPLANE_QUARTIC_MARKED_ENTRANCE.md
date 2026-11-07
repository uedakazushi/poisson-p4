# 平面 H+既約quartic 全源：τ12からの marked-curve 必要条件

2026-11-07、PRIMARY E。producer `exceptional_coverage`、独立確認前。
本稿は type6 nil22-image、base括弧f=0、実五次 D=H*F、Fが既約quartic、HがFの因子でない
全支持に対する幾何学的必要条件を証明する。
全quarticのprojective orbit分類、全marked lineの正規形、既知成分帰属はまだ結論しない。
**F可約の枝へ磁気消去やSaito条件を無断で拡張しない。**

## 1. 元の全sourceからτ(D)=12

`HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND.md` の全24係数f0源を使う。
全二anchor U,VとEuler Eについて det(E,U,V)=unit*H*F であり、
Dは被約なのでSaito criterionによりexponents(2,2)のfree quinticとなる。
従って全sourceで

```
tau(D)=(5-1)^2-2*2=12.
```

ここで本質三変数の平面商を用いること、全磁気六係数が残っていることは元の全源結果への依存。
この必要条件は任意のH,FだけからPoisson sourceの存在を主張する逆向きではない。

## 2. Hとの交点は高々二つ、Fはrational

rを集合としての相異なるH∩Fの点の数とする。Bezoutより交点重複度の総和は四。
各交点でMilnorの積公式を用い、Fの他の特異点も加えると

```
mu(D)=mu(F)+8-r.
```

integral quarticの算術種数は三なので

```
delta(F)<=3,
mu(F)=2*delta(F)-sum_p(branches_p(F)-1)<=6.
```

tau<=muとtau(D)=12よりr<=2かつmu(F)>=4+r。
従ってdelta(F)=3、正規化の種数は零であり、mu(F)は五又は六だけである。
r1ではH∩Fは重複度partition(4)、r2では(3,1)又は(2,2)。
交点はFのsmooth点でもsingular点でもよく、どちらも保持する。

## 3. 漏れのない局所profile九型

delta<=3の平面singularityの枝は次のように数えられる。
unibranchなら、multiplicity二のPuiseuxは(2,3),(2,5),(2,7)、
multiplicity三では最小delta三の(3,4)だけ。multiplicity四以上はdelta>=6で不可能。
これらは順にA2,A4,A6,E6である。
二枝で局所delta<=3なら、二本のsmooth枝の接触I=1,2,3（A1,A3,A5）、
又はA2 cuspとそのtransverse smooth枝のI=2（D5）である。
両枝singularならI>=4に各枝のdeltaが加わり不可能。
三枝はdelta>=3で、等号は三本のtransverse smooth枝D4だけだがmu四なので今回のmu五以上には入らない。

sum(branches_p−1)=6−mu(F)を使うと、全global profileは次の九型となる。

| mu(F) | 全profile | 全枝の条件 |
|---:|---|---|
| 6 | A6 | 一つの(2,7)cusp |
| 6 | E6 | 一つの(3,4)cusp |
| 6 | A4+A2 | 二つのunibranch cusp |
| 6 | 3A2 | 三つのunibranch cusp |
| 5 | A5 | 二smooth枝の接触三 |
| 5 | D5 | A2 cuspとtransverse smooth枝 |
| 5 | A4+A1 | cuspとnode |
| 5 | A3+A2 | tacnodeとcusp |
| 5 | 2A2+A1 | 二cuspとnode |

これはanalytic ADE名に頼ってnodeを除いた表ではなく、deltaと全枝の交点数からの列挙である。
mu五では全体で二枝点がちょうど一つ、他はunibranchである。
全九型の各特異点はsimpleで、tau(F)=mu(F)となる。
global quarticごとのprojective moduliやHの配置は、この局所表だけで有限orbitへ縮めていない。

## 4. rごとの正確なTjurina条件

### r=2

mu(F)=6、mu(D)=12であり、tau(D)=mu(D)。
従ってSaitoのmu=tau criterionにより、**Dの全singularity**がquasihomogeneousでなければならない。
Fは上の四cuspidal型である。
Hの接触partition(3,1)/(2,2)を全て保持し、Hがcuspを通る場合、smooth flex/bitangentの場合を
一つの代表で済ませない。
特にE6 cuspを通る任意のtransverse lineが全てquasihomogeneous pairを与えるとは仮定しない。
この条件は全marked line parameter上で確認すべき制約である。

### r=1,mu(F)=5

mu(D)=12なので再びDの全singularityがquasihomogeneousである。
HがFのsmooth点でintersection四のhyperflexとなる場合を残す。
またHがsingular点を通る場合、四のintersectionを一枝だけでなく全枝へ分配する。
nodeで一枝と三重接触し他枝とtransverseとなる(3,1)、D5でcusp tangentとなる(3,1)、
tacnodeの二枝との(2,2)を、名称だけで同一視しない。

A5の共通接線Hは、Bezoutにより各smooth枝との接触が二、二でなければならない。
Fの二枝相互の接触は三であり、Dの三smooth枝のpairwise接触は(3,2,2)。
weighted homogeneousな三smooth枝のplane germの接触は、適当な重みで全て同じk、
又はtransverse axis一つと他二枝の(1,1,k)しかない。
従って(3,2,2)はquasihomogeneousではなく、この**A5共通接線**caseは除外できる。
A5 quarticと別のsmooth hyperflex lineのcaseまでは除いていない。

### r=1,mu(F)=6

mu(D)=13であり、tau(D)=12だから全Milnor–Tjurina差は正確に一。
F自体の特異点は全てquasihomogeneousなので、差は新しいHとの交点に全て集中する。
Hがsmooth点のhyperflexなら局所DはA7、他の全点もsimpleなのでtau(D)=mu(D)=13となり矛盾。
従ってこの枝でHは必ずFのsingular点を通る。
これはr1のmu六を捨てる理由ではなく、marked intersectionの入口である。

さらにこの枝のHは、A6又はE6の唯一のcusp tangent、又はA4+A2のA4 tangentだけである。
実際、A2のtangent-line contactは三なので、r1で必要な四にならない。
A4/A6では局所x=t²とした接線座標の最初の奇数項が五/七である。
接線の三次項があればA2となるので、接線contactは少なくとも四。
quarticへのBezout上界四により正確に四である。
E6の(3,4)branchでも接線contactは正確に四。
従って3A2の場合はsingular点でもsmooth点でも不可能である。
この段階では三marked型をGa9の全源へ実座標変更したとはまだ主張しない。

## 5. 次に必要な有限whole-family課題

以下はまだ証明していない。

1. 全九profileのintegral quarticを、全projective moduliと特殊orbitを保って正規化すること。
2. r1/r2の各contact partitionについて、全marked Hを分類しtau12及び二次Saito basisを確認すること。
3. 残った全familyの全六磁気sourceを、実可逆変換又はactual Poisson curveで既知成分へ認識すること。

Ga9 power1/power3/linear-wや既知A3+A2 quartic24は有用な入力だが、
この数値表の全marked pairを自動的にそれらへ移したとはしない。
primary arXiv:2407.05819v2 §4.1のQ7例
`x0*(x0²*x1*x2+x0⁴+x1⁴)` は追加のD5+H候補の手掛かりに過ぎず、
全Saito/τ12 sourceの一batch検証前に新しい受理候補や全family証明へ昇格しない。

上の四次marked geometryを揃えるまでは、小weight scan、追加jet、個々の点probeを始めない。
本稿の証明済み結論はr<=2、rationality、全九local profileと記したtau/接触制約までである。
