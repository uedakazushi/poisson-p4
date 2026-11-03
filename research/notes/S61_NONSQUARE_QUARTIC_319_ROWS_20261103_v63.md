# 全重み四の実 row span を319行へ保つ縮約

2026-11-03。開 U0={aFBC!=0}、七 light と一 heavy の正確な weighted normal chart を使う。
d!=0、J!=0、product41、full cubic57 の各点下界はこの row span 等号の仮定に要らない。
未知 core 因子は逆にしない。全四次全射を証明した packet ではない。

## 1. 元の重み四354行

七 light の weight1 と heavy v5 の weight2 からなる全重み四単項式は239個である。
元135式の implicit normal chart では、非選択の weight2 行が6本、weight3 行が21本、
weight4 行が33本ある。これらを適切な重みの単項式で掛ける full Macaulay 表示は

`6*29 + 21*7 + 33 = 354`

行、239列である。33本の weight4 source row は選択34行を除いた

```
12,14,17,21,22,26,28,35,38,41,42,46,48,55,58,63,67,
77,78,88,89,91,96,98,101,109,112,118,119,121,126,128,131
```

である。これら33本はすべて残す。
二次式は公開された q0,...,q5、三次式は全63×91内の正規化 genuine21行を使う。
各正規化が元式の非零有理定数と a,F,B,C の積によることは、既存の独立 row-basis 監査が
全成分で検証している。従ってここで generator を正規化しても U0 の実 ideal は同じである。

## 2. 五つの三次恒等式を七 light で掛ける

既存の自由57行基底の最初の五恒等式は、全三次63行のうち
`48,55,59,57,61` を順に消す。その pivot 係数は正確に `1,1,1,a,−a` であり、
後の恒等式には前に消した行が現れない。恒等式は各々元の全91列で独立受理されている。

各恒等式を各 light vk（7通り）で掛ける。product 行 Ri=qi*vj/ui の寄与は
`(coefficient/ui)*qi*(vj*vk)` であり、174個の quadratic×weight2 行の一つになる。
genuine 行の寄与は、その行×vk である。ui は元の全63行の保存済み単元倍率である。
どの分母も a,F,B,C と非零有理定数だけである。

これで全354行に対する35本の正確な core-rational relation を得る。
消す35行は上の genuine5行×light7である。削除列を同じ順に並べると35×35係数行列は
block 三角形で、対角は各 pivot の七回反復、行列式は

`(1*1*1*a*(−a))^7 = −a^14`。

これは U0 全域の単元である。従って35行を残りの319行で解ける。
残りは元の行の部分集合だから、row span の包含は両方向に成り立つ。

**結論。** 元の全重み四 ideal 部分は、

`6 quadrics *29 +16 genuine cubics *7 +33 genuine quartics =319`

行、239列の表示で同じものを与える。五つの三次関係が特殊化後に rank57 を与えない点でも、
pivot が単元なのでこの行加群の等号は保たれる。

## 3. 保存する有限 circuit と限界

code は全239 column monomial、354 source label、35 relation、35 deleted row、
319 retained row を具体的に列挙する。product から quadratic×weight2 行への積の写像、
各既知単元分母、三角 pivot を exact に検査する。
既存に独立受理された cubic 全91列恒等式を数学的依存とし、そのコピーを新しい四次恒等式へ
掛けた有限操作である。元135式の重い三次・四次展開は再実行しない。

これは新しい、全根を保持する実四次の generator 削減であるが、319×239 が全射であることや、
どの quartic pivot が全開で単元であることも示さない。full57 から四次全射を自動的に推論しない。
heavy-defect 曲線上の既存428×300 proper proof の scope はそのまま保つ。
四次239全射を得た点で初めて、既存の二次・三次 bound と H11 の全底を保持した自由長65定理へ進める。

入口 `code/s61_nonsquare_quartic_319_rows_20261103_v63.py`、
出力 `reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/`。
GB・新しい体・未知因子の局所化・パラメータ行列展開は行わない。
