# M12 実際の一係数による非零性の独立検算

2026-11-05。独立 SOURCE ゲート
`M12_AXIS_ONE_COEFFICIENT_SOURCE_20261105_v63_REVIEW` 通過後、一回だけの
有理数 source 評価と subset DP を行った。親の native code、Bareiss 実装、
producer module は使用していない。

VALUE report:
`reports/S61_NONSQUARE_PREPOWER_M12_AXIS_ONE_COEFFICIENT_20261105_v63_REVIEW/report.json`
SHA256 `c5ebf88f59f93ab7a4a4158ed10f50c89618fa9d1ffb9f8d3e7fab5ef38d5709`。
対象 PRIMARY report の SHA256 は
`1cd0ddbac43a9657d22b951befef4d5834e823e6cc01880a531482a7ad164602`。

67 原子を元の三変数分子・分母から制限し、214 演算すべてを独立に評価した。
全演算の保存された有理多項式との一致、追加23セルの実際の source adapter、
整数化後144セルと既知行 LCD・(a) の冪の一致を確認した。
既存579係数を入力として再使用し、旧572 source 計算は再実行していない。

次数10までの部分行列式を保持する独立 subset DP は、正確に

\[
[a^{15}]\det M_{12}(a,0,0)=
\frac{460045760701204897047977860413366725648080149283864576}
{1756614506664276123046875}\ne0
\]

を与えた。各層の全 mask・係数は `truncated-subset-DP-states.jsonl` に保存し、
層ごとの SHA256 と上界を証明書に記録した。
実際の整数積は249673回、DP の数値スロット上界の観測最大は30592、
source 段階は25730、実行時間は約0.172秒であった。
出力は約1.78 MBであり、事前の資源ゲート内に収まった。

この一係数の非零性により、元の有理関数行列式 (D_{12}) は零関数ではない。
これは29個すべての係数や完全な行列式展開の独立検証を意味しない。
数値パラメータ点・新しい素数・rank search は用いていない。

既に独立検証済みの DIRECT source 恒等式

\[
A\,D_{15}^{\mathrm{norm}}=-T\,D_{12}
\]

と合わせると、旧 chart

\[
\Omega=\mathrm{FULLU}\cap D(\Delta_{407}S_{aux}TN_6\overline N_7)
\]

上のさらに (D(D_{12})) において実際の full57 の十分条件が成立する。
ここで (A) はもともと unit と仮定せず、右辺が unit になることから

\(A,D_{15}^{\mathrm{norm}}\) の両方が unit になる。
DIRECT SOURCE peer は
`S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW/report.json`
SHA256 `f946ca6312c0b5167aa5cc257434d9338727a0ff3a6142496bc2ffa55d3890c8`。
P_K92・h22 の追加逆数や新たな source は不要である。
旧 Ω の補集合を削除せず、全域の full57・quartic239・exhaustion は主張しない。
