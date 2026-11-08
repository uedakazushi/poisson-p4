# 二超平面と既約三次因子：非平面全支持の被覆と残三配置

2026-11-08、PRIMARY `E_reduced_nonnormal`。
本稿は `D=H1*H2*F3`、H1,H2 相異なり、F3 既約の全 rank 四構造を扱う。
**本質三変数の平面退化式を持たない全支持**を MPS40 と既存一つの22次元成分に覆い、
型1+1+3全体の残余を既知の平面十三配置のうち三つだけへ集約する。
型全体の完了、全 E 完了、新成分の追加を主張しない。

新しい三陽性 route は全て別作者独立 PASS。
`TWO_HYPERPLANES_J21_H_X_SS_SEPARATE_REVIEW.md`、
`TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_SEPARATE_REVIEW.md`、
`TWO_HYPERPLANES_P0_T06_SEPARATE_REVIEW.md` に全source/全逆/全曲線/全境界の確認を分離保存した。
本稿の全正常形表と二つの n6零補集合の手帰着は以下に明示し、統合担当の最終読解と区別する。

## 1. 全係数と有限正常形の入口

`TWO_HYPERPLANES_CUBIC_SUPPORT.md` §§1–5 の必要十分入口を用いる。
H1=h,H2=x とし、h-Casimir の全 quadratic 代表を取り

```
P=x∂x∧V+W,  V=(B*y+c*x+d)·∂y,
F3=i_W(dy1∧dy2∧dy3)(V),  actual D=h*x*F3.
```

W の低次項を含む全45係数、元の全 affine Jacobi 四式、逆の projective 代表変更を残す。
`V(F3)=tr(B)*F3` も全係数で成り立つ。
B の半単純部分が非零の全支持は同ノート §§2–4 の全source定理により MPS40 に入る。
可逆/特異、全rank、全binary-form rank、zero-eigen translation の全補集合を含む。

B 冪零は centralizer による全十一正常形へ帰着する。
これは初めに定数項を零に置く仮定ではなく、全 c,d の有限 orbit 分類である。
元の線形45係数解と全残 Jacobi は
`CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json` に固定される。

## 2. 十一冪零型を落とさない全表

座標は `(h,x,u,v,w)`。`S=h+x` とする。

| exact cofactor V | 全 source の結論 | 数学の入力 |
|---|---|---|
| `v∂u+w∂v`, translation L=0 | F3 は w²で割れ、対象なし | 元全線形源と §§5.1 |
| `v∂u+w∂v+h∂w` | 全残 Jacobi から F3 が h または x で割れ、対象なし | 同上の全 n9=0/nonzero 分岐 |
| `v∂u+w∂v+S∂w` | 第一 Jacobi の全線形系が不整合 | 同上 |
| `v∂u` | F3 は v²で割れ、対象なし | 元全線形源 |
| `v∂u+h∂v` | F3 は h で割れ、対象なし | 同上 |
| `v∂u+S∂v` | F3 は S で割れ、対象なし | 同上 |
| `v∂u+h∂w` | 全既約 F3 は F3(h,x,v) の平面入口へ | §3 下記の全 n6 分岐 |
| `v∂u+S∂w` | n6非零の全六源は MPS40、n6零は平面入口 | modular actual curve と §3 |
| `v∂u+h∂v+x∂w` | 全六源は MPS40 | semisimple(3,1,-1),j0 actual curve |
| `v∂u+S∂v+h∂w` | 全三離散源 p=0,±1。p0全五源は T06、p±1全五源は既22成分 | §4 下記 |
| `x∂u+h∂v` | 非平面の全八源/全七源は MPS40、残枝は平面入口 | 既全二 source orbit-closure と二 actual SS curves |

B=0 の translation rank 一では F3 がその共通一次因子で割れ、rank 零では rank 四でない。
最後の rank 二型の H/x swap は八源の同じ全chartへ移り、落とさない。

## 3. 二つの n6零補集合は既に平面入口

`V=v∂u+h∂w` の元19係数源で、残 Jacobi に

```
n6*n16=n6*n17=n6*n18=0
```

がある。n6 非零なら n16=n17=n18=0 なので元の表示 F3 を h が割り、既約対象には入らない。
n6 零なら `F3_u=-n6*h²=0`、さらに全恒等式 V(F3)=0 より h*F3_w=0。
整域上だから F3_w=0、従って F3=F3(h,x,v) である。

`V=v∂u+S∂w` も元の全sourceで `F3_u=-n6*S²`。
n6零なら V(F3)=0 より S*F3_w=0、従って F3=F3(h,x,v)。
n6非零 chart は `(u,v)->(u/n6,v/n6)` の実 PGL で n6一へ正規化でき、
全六自由係数を保持する modular confluence が覆う。
両零補集合を空とせず、同じ既平面13源へ明示的に戻した。

最後の B=0 rank二型でも、F3 の `(hu-xv)` 係数は Ah-Bx。
全 A,B 非零/片方非零は既二 source theorem が覆う。
A=B=0 は F3(h,x,w) となり、同じ平面入口へ戻る。

## 4. 三つの actual whole-source 橋

`TWO_HYPERPLANES_J21_H_X_SS_CONFLUENCE.md` と別作者 review は、
全六源を source 分母の無い三角 PGL で二 parameter (s,k) へ送り、
半単純 `(3,1,-1)`,j=0 の全A3親から全 s,k の literal curve を与える。
親の全像を含む ambient 成分が SNC 点に会うので全親は MPS40 に入る。
元の全150 normalizer/逆、全曲線350、Pfaffian と全端点が両作者で一致する。

`TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_CONFLUENCE.md` と別作者 review は、
n6一全六源を source 分母の無い PGL で k=n16 全モデルへ送る。
actual polynomial `Qe+k*Ve∧Ze` は generic e で cofactor diag(0,1,1) と零 kernel translation を持つ。
§4の既受理全非零半単純支持定理と閉性により全 k endpoint が MPS40 に入る。
必要な generic chart の分母は e と 1+2ke だけで、端点は全 k に正則。

`TWO_HYPERPLANES_P0_T06_CONFLUENCE.md` は p0 の元五自由係数を全逆三角 PGL で l=n3 全モデルへ送り、
真 T06 source に

```
p0=e+2e²l, p1=-e-4e²l, t0=e, t1=-e
```

と actual cubic-preserving moving marked-pair PGL を代入する。
Euler 商全150係数は e=0 に正則、端点は全 l の元モデルに一致し、全350も零。
全三 route で parameter/磁気係数零は同じ literal family に保持する。
全三 route の独立 review は完了しており、全 source の PGL 飽和閉包まで MPS40 に入る。

p=±1 二源は既受理 `R4_TWO_HYPERPLANES_CUBIC_22` の同じ全22次元成分を覆う。
両指定点で全D103/実source47により smooth unique22、実 PGL/Euler 同定で同一成分と証明済み。
今回その旧成分を重複加算せず、旧 source/独立 Fraction component review をそのまま依存にする。

## 5. 型1+1+3の正確な残余は平面三配置

上表で平面に戻る枝は actual D が本質三変数である。
既受理 essential-three全cotangent/Saito theorem と
`TWO_LINES_PLANE_CUBIC_SUPPORT.md` の full13 finite-source theorem により、
全磁気六係数を含む十三配置へ尽くされる。
うち九源は既23次元成分、節点三次二節点接線は全六磁気の toric collision で MPS40 に覆われる。
従って未帰属として残すのは正確に

1. `node_flex_tangent_through_flex`：節点三次、変曲接線と変曲点を通る通常接線。
2. `cusp_flex_secant_through_flex`：尖点三次、変曲接線と変曲点を通る一般割線。
3. `cusp_cusp_tangent_cusp_secant`：尖点三次、尖点接線と一般尖点割線。

各 full6 source をそのまま残す。旧二次/固定三次の零、ある部分的 parent 陰性から未帰属を消さない。
本稿で追加する成分は零。証明された全source被覆と、残三源の actual membership 義務を区別する。
