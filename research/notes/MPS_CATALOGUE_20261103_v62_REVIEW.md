# MPS の一般有理族と支配性の独立点検

2026年11月3日。対象は v62 の raw catalogue、有限 Maurer–Cartan 補正、
PGL・period・edge 微分の三証明書である。
**T02–T38 の37個の一般有理 Poisson 族と、T01–T38 の対応する像次元を受理する。**
既知成分のラベルと成分性には下記の MPS 局所定理と既存の図式照合を依存として明示する。
全51階数四成分の S61 定義 ideal が完成した、あるいは大域尽くしが示されたという意味ではない。

## 一般係数と独立な検算

4×4 biresidue 行列 B の `−Pf(B) B⁻¹` は、実際には B の成分に線形な Pfaffian 双対である。
Pf(B)≠0 の領域で、Poisson 係数と二形式の period 係数を混同しない変換である。
装飾辺の monomial を加えた raw 族19件は直接 Jacobi を満たす。
残り18件の失敗を保持し、全 period 有理関数体上の有限補正を別に構成している点を確認した。

補正の列は toric Jacobi 微分である。`Jac(toric)=0` なので
`Jac(toric+V)−Jac(V)` はちょうどその線形微分になる。
辺重みが独立な場合、150個の二次双ベクトル monomial の各重みには高々一つの
非負 smoothing 多重次数しか対応せず、候補が有限である。
ただし候補の列挙を尽くし定理とはせず、最後の全 projective Jacobi の消滅を結論の根拠にする。
周期のある T21 は raw の恒等式で処理され、独立重みの仮定を誤用していない。

監査用の `code/v62_independent_review_20261103.py` は作者の二つの Jacobi 関数を呼ばず、
保存された最終係数だけから有理関数体上の sparse 多項式環を作る。
20個の三添字 Schouten 成分を直接微分し、Euler ベクトルと wedge した全四添字成分を確認した。
T02–T38 の37族すべてで厳密に零となった。
さらに各最終項の空間 torus 重みが smoothing 多重次数の辺重み和と一致することと、
period・smoothing の同時スカラーに関する次数1の斉次性も独立に確認した。

## 像次元と非零小行列

作者が保存した38個の有理微分行列について、指定部分行列の係数を分子・分母から
F32003へ正しく還元し、行交換・除算を別実装して全38行列式を再計算した。
全て保存値と一致して非零である。
これは QQ 上の階数**下界**であり、標数零の階数上界を有限体から推測する操作ではない。

18個の補正族については、さらに全最終係数の分母がこの**同じ保存period標本**・辺変数0で
非零となることを独立に検査した。完成族とraw族の差の値、および全period・辺変数による
一次微分がその点で零になることも直接確認した。従ってraw族から作られた保存微分行列は、
実際の完成Poisson族の微分行列でもある。証拠と全分母値は
`reports/MPS_CORRECTION_DENOMINATORS_20261103_v62_REVIEW.json`、実装は
`code/mps_correction_denominator_review_20261103_v62.py` にある。

上界は直接にも得られる。確認済みの torus 重み共変性により、対角 T4 は
`PGL5×(period,edge)` の源に自由に作用し、像を固定する。
従ってアフィン像次元は高々 `24+s+e−4=20+s+e`。
上の非零 minor がこの下界を達成するので、像次元は正確に `20+s+e` である。
スカラー斉次性から像閉包は錐で、射影次元は `19+s+e`。
S61 との既存の次元差9から、対応する S61 成分次元は `10+s+e` になる。

## 成分ラベルに使う一次文献の仮定

MPS, *A local Torelli theorem for log symplectic manifolds*,
[arXiv:2010.08692v2](https://arxiv.org/html/2010.08692v2) の一次本文を今回再確認した。

- [Theorem 6.16](https://arxiv.org/html/2010.08692v2#S6.Thmtheorem16) は、
  smoothing subdiagram ごとに対応する局所既約成分の一意性を与える。
- [Corollary 6.23](https://arxiv.org/html/2010.08692v2#S6.Thmtheorem23) は、
  resonant edge の所定の H1 消滅などの条件で deformation fiber が滑らかであることを述べる。
- [Example 6.24](https://arxiv.org/html/2010.08692v2#S6.Thmtheorem24) は、
  座標超平面境界を持つ射影空間では同 corollary の第二条件が自動で成立すること、
  P4 では Table 2 の40族の大域的相違を直接確認したことを明記する。
  これが [Theorem 1.4](https://arxiv.org/html/2010.08692v2#S1.Thmtheorem4) の40成分に接続する。

既存の `MPS40_CROSSWALK_20261031` は、period 制約、装飾辺とその次数、
基点での全 smoothable edge のリストを Table 2 に照合している。
新族は同じ一般 period 源を持ち、各対応 edge の非零な一次 cocycle を持つ。
上の局所定理により該当する局所成分へ入り、像次元が一致するため既知成分を支配する。
単に次元表が同じという理由だけで内部IDを同定してはいない。

文献の外部証明全体や40 SVGの照合そのものを今回別実装したとはしない。
一般 S61 rational chart の正規化と134式の検算は、作者の後続の有限 chart 証明書を別途点検する。
今回の監査記録は `reports/MPS_INDEPENDENT_20261103_v62_REVIEW/report.json`。
