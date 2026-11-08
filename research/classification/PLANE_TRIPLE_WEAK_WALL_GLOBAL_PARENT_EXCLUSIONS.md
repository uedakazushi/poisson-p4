# 二つの E1DX weak wall は global T39/T40 の外にある

2026-11-08。PRIMARY `E_reduced_nonnormal`。producer exact arithmetic と既受理の大域評価-fiber 定理による結論。独立算術確認は未実施。

## 命題と全源の範囲

`PLANE_M3_TRIPLE_WEAK_WALLS_TWO_ACTUAL_ORBITS.md` の実六次元 unipotent 群と全逆により、τ=4/21 および τ=1/12 の各 complete affine W六源はそれぞれ一つの actual PGL5 軌道である。本稿はこの二軌道の **全点が global T39 および global T40 のどちらにも属さない** ことを示す。他の二十四次元成分の境界への所属、新成分性、全分類の完了は主張しない。

## 固定された一つの rank 四点と実 frame

元座標を(a,b,c,z,w)とし、上記 note の canonical particular tensor を用いる。両点とも p=[0:0:0:0:1] で actual divisor は非零定数 E を持つ。E=1/2 (τ4/21)、E=−1 (τ1/12)。このため p は各軌道の構造から最初に固定できる symplectic point であり、point/frame 探索は要らない。

S61 座標への actual PGL frame はそれぞれ

```
(x0,x1,x2,x3,x4)=(w,a,b,c,2z),
(x0,x1,x2,x3,x4)=(w,a,b,c,−z).
```

p の affine spatial skew matrix はこの frame で J12=J34=1 に一致する。両行列と真の逆による full150 coefficients の pushforward 往復を report に保存した。perspective normalization は両点とも恒等である。Euler quotient、全五十五 triangular elimination の再構成を含む全150係数、元 reduced134式、projective Jacobi が `S61.verify_point` で一致した。

## 全四十五 graph row の正確な値

両 canonical tensor の全61座標は report に保存した。非零は零始まりの位置 3,5,28,46,54 のみである。

| τ | u[3] | u[5] | u[28] | u[46] | u[54] |
|---:|---:|---:|---:|---:|---:|
| 4/21 | −7 | −6/7 | 15/7 | 3/14 | 3/7 |
| 1/12 | 20 | −3/2 | 3/2 | −3/4 | 3/4 |

全十六 weight 一座標から計算される T39/Q51 および T40/Q52 の quadratic section を引いた全四十五 residual は、**どちらの graph についても**位置5,28,46,54だけが上表の同じ非零値となり、残四十一 residual は零である。従って両 S61 点は両閉 graph の外にある。

`S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md` と独立 `S61_GLOBAL_COMPONENT_FIBER_SEPARATE_REVIEW.md` は、connected Sp4 評価 stabilizer と global Ga4 perspective product により、全 global T39/T40 の S61 交点が各閉 A16 graph に **全点で**一致すると証明している。generic beta、非零 covariant、normal first jet の仮定はない。このため上の陰性は特定の normalized chart の像だけの除外ではなく、全 global component の除外である。PGL5 不変性と全 W六 transitivity を戻すと命題を得る。

## 証明書と残る幾何

実行前 gate は PASS。新 mode は reusable runner の旧五十 FunctionDef AST を保ち、公開9468048fでの source全文・helper全文・当時の current契約を preappend archive に固定した。完了済み G 契約の binding を新しい E mode へ更新していない。

`CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_GLOBAL_GRAPHS_20261108.json` に二つの full Pi、actual D、固定点の値、frame・全逆、全61座標、全45 residual、full150/134検証を保存した。4.795秒、最大RSS150100KiB。これは新 Jacobi rank、追加 minor、jet、点列の計算ではない。

両 tensor は absolutely irreducible essential 五の actual D と二十三次元 actual orbit を持つ。既二十三次元 nil22long/K4a 軌道閉包とは singular plane 数で異なるが、二十四次元 MPS 親への boundary 包含は未決である。T39/T40 陰性から新成分と結論してはならない。次の最弱候補は同じ二つの全軌道における full ambient tangent の一回だけの判定である。D102 なら二十三次元 source 閉包の全成分性へ、D101 なら適切な非零 weight 法線と完全有限 equivariant slice から二十四次元親への実曲線へ進める可能性がある。まだこの次計算は実行していない。
