# 真の通常三次初期関係：完了した有限次数証明書

[証明と範囲](../../research/notes/Q_ORDINARY_CUBIC_INITIAL_20261103_v63.md)を参照。
最終結果は `summary.json`。抽出直後の `report.json` にあるrank未計算フラグは抽出phaseの状態であり、最終QQ/GF rankは `exact-rank/report.json` と `modular-rank/report.json` に保存している。

- 全350三次jetをQQ上で抽出し、physical Jacobiと照合。
- 57個の二次部分の定数消去で135本の非零な真の三次初期式を得た。
- 二次式の全35変数倍1995行も含めると、QQ上でrank1894→1911、真に新しい三次classは17個。
- QQの通常Hilbert関数は次数3まで `(1,35,573,5859)`。
- 唯一の座標13平面上のthrough3制限イデアルは `(z8^2*(12*z14+z17+6*z3))`。台は二つの12平面であり、全局所環の次元12を意味しない。

`cubic-initials.json.gz` は全QQ係数と全定数消去操作、`cubic-blocks.in.gz` と `integer-blocks.in.gz` はnative行列入力を保存する。後二者を解凍したSHA256は各phase reportにある値と一致する。`plane-certificate.json` は別の短い実装による全350行操作の再生と正確な因子式の検査結果である。

再現コマンド、V=Z_(32003)上の有限次数の意味、未解決の高次関係は上記証明ノートに明記した。大きなGBは実行していない。
