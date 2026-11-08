# 本質三変数のH+二conic：全配置・全磁気源の被覆

2026-11-08、PRIMARY `E_reduced_nonnormal`。
全配置の尽くし、旧五源の被覆、新chord及びtan-simpleの全実衝突は、
それぞれ独立review PASS。T16/T17へのnamed同定も別作者が全係数で確認した。
本稿はこれらの受理済み入力を合わせた全source被覆定理であり、新成分数は加えない。

## 1. 被覆対象と全源入口

一般階数四のprojective Poisson構造で、実退化因子が

```
D=H Q1 Q2
```

と分かれ、Hは単純超平面、Q1,Q2は相異なる既約quadric、Dの本質変数数が正確に3とする。
canonical essential quotientのP2では二quadricはsmooth conicsである。
全source還元 `TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` §§1–4と
独立 `TWO_QUADRICS_PLANE_CONICS_REVIEW` / `PLANE_SOURCES_REVIEW` を使う。

その入口は選んだtensorの最高次部分だけを扱わず、全fiber係数と全Jacobiを残す。
全rank四sourceは二次Saito frame U,Vと任意の六係数二次磁気項Rから

```
P=U∧∂z+V∧∂w+(−az−bw+R)∂z∧∂w, [U,V]=aU+bV
```

の形となる。逆にこの全六磁気sourceは元のPoissonと実Pfaffianを保つ。
frameの変更は可逆なfiber GL2へ吸収されるので、一つのbasisに限った分類ではない。

## 2. 配置と全補集合

Bézoutと局所Milnor budget、free exponents(2,2)のTjurina数12から十候補へ尽くす。
三候補は非零な一次gradient syzygyにより全域でfree(2,2)と両立しない。
残る七配置では、全18二次gradient係数の完全kernelがrank2、
Saito determinantが非零であることを独立に確認している。

唯一の配置parameterを持つshared-tangent一般族では、同一conic又は可約conicの境界を別保持し、
source chartのu=1/2は別の非零minorで同じrank2frame bundleへ含める。
三候補の排除はgeneric pivotの零点に依存せず、全parameterで成り立つ一次syzygyである。
したがって未検証のparameter wall又は未知frameを黙って省いていない。

## 3. 七つの完全源と閉親

| 完全配置 | 全源の被覆先 | 証明の入口 |
|---|---|---|
| shared_tangent_generic | CPLUS_EXT | 全frame bundleと既知唯一成分点、全u/全六磁気 |
| shared_tangent_contact3 | CPLUS_EXT | q2(t)=q+xy+t*xz、全rank2kernel/frameの正則延長 |
| contact31_common_bitangent | R4_PLANE_TWO_CONICS_BITANGENT_23 | 真48切片の非零一法線二次支持と独立全ambient算術 |
| contact4_other_tangent | R4_PLANE_TWO_CONICS_CONTACT4_TANGENT_23 | 真48切片の一法線三次支持と独立全ambient算術 |
| contact4_secant_at_contact | 既知beta0有限閉和 | 全六磁気の正則translation/Darbouxで実J+Qへ |
| contact31_chord | T16 | 全六磁気のactual SNC confluence、別作者PASS |
| contact31_tangent_at_simple_intersection | T17 | coupled actual SNC confluence、全六磁気・別作者PASS |

旧五源の詳細は
`TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md`、
`TWO_QUADRICS_PLANE_COMPONENTS.md`、
`TWO_QUADRICS_PLANE_EVEN_SUPPORT.md`。
新二源の全式は
`TWO_QUADRICS_CHORD_SNC_CONFLUENCE.md`、
`TWO_QUADRICS_TAN_SIMPLE_SNC_CONFLUENCE.md`。
T16/T17のfull completed sourceとのexact Euler/parameter同定は別の
`TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` にあり、元の有限閉和被覆はその番号同定に依存しない。

具体的な独立入力は次の通りである。

* 全十候補・七free source：`CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS`、同`CONICS_REVIEW`、
  `PLANE_SOURCES`、同`PLANE_SOURCES_REVIEW`の各20261107 report。
* CPLUSの全24成分と唯一成分点：`legacy/v08/sections/cplus_component.tex`。
  その点から全frame bundleを覆い、contact3まで実延長する論理は
  `TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` §5。
* 二つの23成分：`TWO_QUADRICS_PLANE_COMPONENTS.md`、
  `TWO_QUADRICS_BITANGENT_SEPARATE_REVIEW.md`、`TWO_QUADRICS_OTHER_TANGENT_SEPARATE_REVIEW.md`。
  全D/Fと真切片・二次/三次の算術は`PLANE_AMBIENT/NORMALS/CUBIC_SUPPORT`と
  `PLANE_COMPONENT_REVIEW/PLANE_CUBIC_REVIEW` reportに保存される。
* secant-at-contact：`TWO_QUADRICS_PLANE_EVEN_SUPPORT.md`と`PLANE_EVEN/PLANE_EVEN_REVIEW`。
  そこからの全beta0被覆は`BETA_ZERO_COVERAGE.md`の既存定理を入力とし、
  同定前のQ45全22台定理は`legacy/v07/quadratic_classification.tex`のTheorem Q22である。
  本completionでその旧Q45定理全体を再証明したとは主張しない。
* 新二実衝突：`TWO_QUADRICS_CHORD_SNC_SEPARATE_REVIEW.md`、
  `TWO_QUADRICS_TAN_SIMPLE_SNC_SEPARATE_REVIEW.md`と各20261108 report。
  後者はnamed T16/T17の全五parameter・全六括弧・Euler gauge・逆bindingも独立に確認した。

## 4. 新二源のactual-family論理

各親は五parameterの既約な実Poisson源であり、全Jacobiを満たし、実SNC rank四点を含む。
親の全像閉包を含む一つのambient成分はそのSNC点も含む。
採用済みのSNCに会う全成分尽くし定理により、その成分はMPS40の一つである。
従って親全体が既知有限閉和へ入る。
SNC一点が既知成分に入るという事実だけを逆向きに外挿していない。

source/periodとその積関係を保つ実PGL行列（detは±epsilon²）とtensor epsilon倍から、
負のepsilon冪が全て消えた多項式族を作った。
rank4の磁気translationと親二periodの商rank2を合わせた六写像は可逆で、
任意の全六磁気係数を端点に一致させる。
全150端点係数、元350/projective Jacobi、strict10、実Pfaffian、全座標及び磁気逆を保存した。
epsilon≠0は既知親の真の点、epsilon=0は全対象sourceなので、閉性が全点所属を与える。

磁気parameterは一つも反転しない。磁気零、任意の零係数、途中のperiod rank低下も同じ実族に含まれる。
一法線二次・三次が非孤立だった旧計算は、新しい正橋の仮定ではない。

## 5. 完了条件と残る範囲

七つの完全源の有限閉和被覆が独立受理まで揃った。
全source入口が尽くすことから、§1の**本質三変数H+二相異既約quadricの全点被覆**が結論する。
新component数を加えず、全PGL飽和及び全磁気境界を含む。

これはH+二quadricの非plane部分、lambda零pure-nil等の別source、
他の退化因子型、又は全Poissonスキームの尽くしを完了したという主張ではない。
T16/T17等の内部でどの全境界が同じ軌道か、非被約schemeの精密構造も要求しない。
