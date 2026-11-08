# 残る二conic平面源：二つの24次元親を除外する

2026-11-07、PRIMARY `E_reduced_nonnormal`。純粋なmarked divisor幾何による補題。
全周囲成分性、新成分、他の親の不存在は主張しない。
全源は `TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` とその独立source reviewの既済入力である。

## 対象とproper factor flag

q=xz−y²とする。残る二源は

```
contact31_tangent_at_simple_intersection: H=z, (q1,q2)=(q,q+xy),
contact31_chord: H=y, (q1,q2)=(q,q+xy).
```

両方とも相異なる二つのsmooth plane conicsと唯一の単純Hを持つ。
六磁気係数は実D=H q1 q2を変えない。
exact essential3開では方向微分kernelから商平面が正則に復元できる。
H及び二conicのprojective factor flagはproperなので、以下の閉marked条件は
元のPGL5成分の閉包内でも保たれる。二conicの交換は許すが、Hの別選択は存在しない。

## Ga9 cuspidal plane-quartic24の除外

`GA9_CUSPIDAL_MARKED_DISCRIMINANT_FRONTIER.md` §§1–2の最も弱い必要条件だけを使う。
同成分の一般marked因子はD=H F4で、F4|Hは一点4pである。
四乗binary formの集合はproper乗法像で閉なので、F4が可約になる境界でも
F4|Hが非零なら四乗条件が必要である。quartic二次coneやdiscriminantは今回は不要である。

tan-simpleではz=0への制限が

```
q1 q2|z=0=−y³(x−y),
```

相異なる二点の重複度3+1であり、四乗ではない。
chordではy=0への制限が

```
q1 q2|y=0=x²z²,
```

相異なる二点の2+2であり、やはり四乗ではない。
従って両全六磁気源は `R4_GA9_CUSPIDAL_PLANE_QUARTIC_24` の外にある。

## CPLUS_EXTの除外

`TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` §5及び既済CPLUS_EXT成分定理によれば、
CPLUS_EXTの一般marked因子はHと二smooth conicsであり、Hは両conicに同一点で接する。
H上で各二次形式がlinear formの平方となる条件は閉である。
同じ接点という追加条件を使わず、より弱いこの二つの平方条件だけで十分である。

tan-simpleではq1|z=0=−y²は平方だが、q2|z=0=y(x−y)は二つの相異なる根を持つ。
chordではq1|y=0=q2|y=0=xzで、両方とも相異なる二根を持つ。
従って二conicを交換しても必要な二つの平方条件が成立せず、両全源はCPLUS_EXTの外にある。
target conicはいずれもsmoothなので、途中のrank落ち境界やfactor chartの選択で
このproperな閉接線条件を免れることはできない。

## 適用範囲

これは磁気一点の陰性でも、一般点のincidenceから特殊点を外挿した議論でもない。
同じ実Dを持つ全六磁気係数、磁気零、全PGL飽和に一様に成り立つ。
MPS T16/T17/T36/T37/T38等の別の親は残り、全E被覆は未完である。
旧二法線二次・三次の非閉鎖計算を再実行していない。
