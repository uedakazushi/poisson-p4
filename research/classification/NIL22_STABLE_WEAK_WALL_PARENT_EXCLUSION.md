# Stable 二全源は新 weak-wall 二十三成分へ全体として入らない

2026-11-08、PRIMARY E、既済入力からの純手帰結。新計算はしない。
個々の磁気特殊点の非所属ではなく、whole-source parent 候補の必要な除外だけを主張する。

stable_line_conic と stable_three_lines の全六磁気源の PGL5 飽和像閉包を Z とする。
source は affine 六空間なので Z は既約。
既存全 source 行列は両方で rank48、Euler24 と非零 scale 一を引いて
`dim Z >=23` を与える。
全 actual divisor は三座標だけで書けるため、Z の全体で essential variables 数は三以下。
これは五つの degree四 partial derivatives の coefficient matrix の rank三以下という
閉条件であり、PGL5 作用と closure に保たれる。

新二 weak-wall 成分 C は、それぞれ全 ambient component として dimension23。
その dense actual orbit の divisor は絶対既約かつ essential-five。
もし Z が C に含まれるなら、既約閉集合の dimension から Z=C。
すると C 全体が上の essential≤3 閉条件を満たし、dense essential-five divisor と矛盾する。

従って **どちらの weak-wall23 成分も、stable 二全源のいずれをも覆わない**。
各 source 上の C の逆像は proper closed なので、その非所属開集合も非空。
特殊磁気 subsource や個別点が C に入る可能性は残す。
この結果だけで別の二十四次元親の存在、stable の成分性、又は全点非所属は推論しない。

数学入力：

- `NIL22_STABLE_ORBIT_FRONTIER.md` §1 と元 `CLASSIFICATION_NIL22_IMAGE_STABLE_CUBIC_20261107.json` の全作用/source rank。
- `PLANE_TRIPLE_WEAK_WALL_FULL_COMPONENTS.md` の dimension23 と dense essential-five orbit。
- `PLANE_TRIPLE_WEAK_WALL_COMPONENTS_SEPARATE_REVIEW.md` の独立受理。
