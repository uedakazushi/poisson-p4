# v62 非平方 core：二回目の探索引継ぎ

2026-11-03。旧 bundle は保持し、追加の41ファイルを `nonsquare_core_wip_v2.tar.gz` とその manifest に保存した。これらは探索記録である。証明済みの新結果は別の `S61_BINARY_NONSQUARE_CUBIC_REPLAY`、`HEAVY_CURVE_CUBIC`、`HEAVY_CURVE_QUARTIC` の v62 report/note を参照する。heavy 欠損曲線の周囲除去はそこで証明したが、全 marked 非平方層はまだ尽くしていない。

保存した未完了計算は次の通り。

- 六二次式の light 一次倍：既知の一つの syzygy を除いた41×84行列は F32003 では全 aFBC 開で階数41。QQ の module Gröbner 計算は600秒で時間切れだった。既知単元21 pivot 後の20×63表示も240秒で時間切れ。QQ の全点定理は未証明であり、WIP report は `all_passed=false` とした。`marked_square` agent が定数行・列結合で隠れた単元を探す別表示を継続している。
- B-free の四次 prolongation：元135式から source を再構成し、正規化 rows の SHA が `32bac22b558bf2267506dd34414ba76ab3b5382acddd01c5f38b6948f089bfc1` と一致した。321×238行列は約44MB。素朴な既知単元 pivot は32段で係数項数上限に達した。Markowitz 型の非零数による候補選択で51段まで進み、270×187・約476万 QQ 係数項の残余となった。この大きさを小さい annihilator の定理へはまだ変換していない。
- a=1 の直接 t=1 模型は、G6 を含めても F4 の有限体計算が6GB virtual-memory cap で `std::bad_alloc`。未知因子を除かず既知 a,d,F,B,C content だけを除く前処理でも改善しなかった。同じ入力の無制限再実行はしていない。

大きい再生成可能な B-free 行列・Schur file は bundle に重複収録せず、manifest にサイズと hash を残した。`fresh-bfree.py`、`bfree-unit-cover.py`、native Schur source とログから再生成できる。必要な native compile は公開済み Cloud の `c++`、`pkg-config flint`、`-lgmp` を使う。元のセットアップをやり直す必要はない。

Macaulay2 1.24.11 の互換性に関する追加所見を保存した。`gb(..., Algorithm=>LinearAlgebra)` は小さい QQ 2×2 module でも `array out of order` の後 SIGSEGV になった（`la-test.m2`, `la-test.log`）。本研究の全証明書は既定の Inhomogeneous `gb` を使っている。また、この版の `Core/gb.m2` では `groebnerBasis(..., Strategy=>"F4")` は prime finite field のときだけ MGB/F4 engine を呼び、QQ では通常の `gb` に戻る。従って QQ の実行を F4 実行と表示してはいけない。重み付き有限長を測るときは `degree` とベクトル空間長を混同せず、`numColumns basis` を使う。

次は、heavy 階数7の場所に残る genuine cubic/quartic 欠損曲線を具体的な小さい行列または proper な有限証明書で覆う。今回の heavy 曲線と同様、有限体の点を見つけるだけで QQ の全曲線除去とはしない。
