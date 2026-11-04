# H15 有限候補上の実 quartic239：独立レビュー

2026-11-04。**FINAL PASS**。
対象 primary は `S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63`、
report SHA-256
`7041b20efce46e4a4f2596cc64a7f313bb37f3ae66be8943dcc857e061976291`。

結論は、許容 P100 曲線上の **H15 の全有限候補**で、実319×239 quartic
係数行列が239列を生成することである。H15 の外での239、P100 曲線全体の
239、全 FULLU の239をこの receipt からは主張しない。

## 独立に確認した固定入力

著者の coefficient/NF evaluator や rank search は再実行していない。
別実装の token checker によって、`config.txt`, `shared.txt`, `cubics.txt`
を元の凍結 QQ 証明書から順に再構成して全 token を比較した。

- 原 retained319 catalogue の quadratic174/cubic112/raw quartic33を全保持。
- 完全16 cubic×91係数を圧縮57×91から逐次照合。
  148858 core terms と decompressed source SHA が一致。
- Oracle 全33 raw 根、Tden による変数置換、Qhat の DH^-2、座標逆変換49。
- NF122 の全 scalar/linear closure、H² pair77/131の u,v,U。
- actual delta7 frame の112×23右逆 L の全 active42行。
- 176101個の実 QQ 係数について分母の32003-unit性を再検査。
  別に H の16係数も検査した。

この token 照合は新しい parameter polynomial の展開や source query を
含まない。原係数を有限体入力へ移す全写像を検証している。

## source frame と全 pure tail の監査

実装を source catalogue と受理済み FrameAdapter/NF/H1右逆の恒等式に照らし
監査した。物理 Tn 多項式を joint 座標に移してから NF117 を適用している。
原16 cubic の112 products は元の正規化 frame のまま保存し、Section を
作る部分だけに actual delta7 源変換を入れている。

\[
\Gamma=f_0-(f_2/U_2)\chi_0
 -(f_1-(f_2/U_2)\chi_1)
   (L_{\mathrm{act}}^tD_{\delta7}\Pi_{\mathrm{orig}}),
\quad \delta=-64000C^2,\quad U_2=(50aBC)^2.
\]

\(\chi\) は raw77/131 から構成され、各因子で H²係数=Uを検査する。
23個の section identity は受理済みの whole-curve 右逆の特殊化であり、
今回の著者 native でも全因子で実測検査してから L/cache を解放している。
raw \(f_0\)、\(\chi_0\)、全 cubic pure tail は保持される。
21個の pure quartic roots や未構成の70+15+4核の完全性は仮定しない。

## 独立の有限体・固定 rank 証明

係数配列の別実装で次を確認した。

1. Hbar の全因子積と H の QQ 係数からの直接 reduction が一致。
   次数 `4,1,1,1,2,3,3` の相異なる monic 因子を全保持。
   Frobenius chains と四個の必要な Rabin-Bézout tests を再検証。
2. 五つの既知曲線 localizer が全因子で非零。
3. native の binary log の **全21329逆元積**を独立乗算し、積1を確認。
   source denominator と rank pivot の双方を含む。
4. 全七 Γ binary は元319 source IDs を過不足なく保持。
   各 selected98×98行列の全係数が対応する実 Γ binary 行と一致。
5. 著者の固定 pivot 列だけを使い、別の C++ field-array 演算で全98段を検査。
   pivot search、rank search、別の行/列選択はしていない。

原列順の determinant は全七個とも著者の保存値に一致した。

| 因子次数 | 独立再計算した非零 determinant の係数列 |
| --- | --- |
| 4 | `[8688,3754,1861,11981]` |
| 1 | `[31268]` |
| 1 | `[27666]` |
| 1 | `[6825]` |
| 2 | `[17263,9304]` |
| 3 | `[326,21025,20839]` |
| 3 | `[19211,2000,22829]` |

## finite integral lift と結論の範囲

\(V=\mathbb Z_{(32003)}\), \(S=V[z]/(H)\) は monic H により有限自由15。
必要な源係数・frame・section 分母の全因子上での逆元が検証されているので、
有限 Nakayama により、それらは S の単元である。従って実 Γ は S 上に
定義され、全特殊ファイバーの CRT と七個の固定非零98-minor により
\(\bar\Gamma\) は全射となる。cokernel は有限 V-module で mod32003 が零
なので、Nakayama により零。標数零への移行で体性・既約性・被約性を
仮定したり、因子や点を捨てたりしていない。

受理済み117個の quadratic unit pivots、実 H² unit 一個、実 H1 section23
を戻せば、\(117+1+23+98=239\)。これが **全 H15 候補上**の実 source-bound
quartic239を証明する。H15 外の全曲線239は別問題である。

## 再現・固定情報

review code:
`code/s61_nonsquare_P100_H15_full239_finite_review_20261104_v63.py` と
`code/s61_nonsquare_P100_H15_fixed_rank_review_helper_20261104_v63.cpp`。

review report:
`reports/S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63_REVIEW/report.json`
SHA-256
`1125cf14c2594385395aee6a397e66acd9616230a3c4989ce709fcfdf6df1076`。

全2.333秒、独立 fixed-rank native0.137秒。
出力は別 review directory のみ。著者の凍結入力・コード・出力は不変。
