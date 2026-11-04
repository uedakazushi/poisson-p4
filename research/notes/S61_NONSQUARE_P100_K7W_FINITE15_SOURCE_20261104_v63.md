# K7W を第九行とする有限 H15 検査の入力固定

2026-11-04。係数計算前に固定した source-only 記録である。

旧五個の 9 次小行列の正規化後の厳密な共通因子を、既存証明書の
monic な (H\in\mathbb Q[z]), \(\deg H=15\) とする。今回の入力は
旧 cover case9 の八行 `[9,2,3,5,7,6,1,4]` に、**新しい**原319行列の
行226、正規化 cubic C68 と physical \(v_3\) の積を加える。
physical \(v_3=W\) は親 K7W 記録で検証済みである。

九列は old local98 の `[0,1,2,3,4,5,6,9,28]`、original122 では
`[0,1,2,4,5,6,7,10,33]` である。全九列を含めた各行の LCD は、
P100 曲線の既知の五つの unit
\(z,z+20,z+60,203z+2700,z+12\) の冪だけから作る。
入力は1057 RF terms。全 QQ 分母と monic H の係数は32003-integral。

有限環は \(V=\mathbb Z_{(32003)}\), \(S=V[z]/(H)\) であり、V 上自由15。
体・既約性・被約性・平坦性の追加仮定は置かない。検査対象は
**全特殊ファイバー**であり、曲線の localizer による点の除外をしない。

この receipt では determinant、remainder、rank を計算していない。
旧 K2W の失敗した determinant やその入力データは読み込まず、後続計算で
再利用するのは汎用の exact Bareiss エンジンのコードだけである。

固定された出力:

- `reports/S61_NONSQUARE_P100_K7W_FINITE15_SOURCE_20261104_v63/report.json`
- `source.json` SHA-256
  `e5d2efd5516df1da4dfc90e67220671f2cd4dfea1ab511e58382087d05dfd9c1`
- `source-input.txt` SHA-256
  `7cb73e448b92ff76646be5d725c7b400a843942e94aee06086070e3a529c4814`

コードは `code/s61_nonsquare_P100_K7W_finite15_source_20261104_v63.py`。
この入力を用いた一回の計算は、別記
`S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63` に保存する。
