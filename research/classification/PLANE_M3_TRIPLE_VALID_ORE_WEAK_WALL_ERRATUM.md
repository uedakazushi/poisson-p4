# corrected triple Ore：二つの consistent wall の弱い G 結論

2026-11-08、独立 reviewer quartic_sources の読解により確認。
immutable 十八 module report の数値・行列・kernel・G flag は正しい。
本訂正は producer aggregate の強すぎた**空性**の文言に対するもので、
G 非優越の結論は変わらない。
旧 report／凍結 proof／manifest は書き換えない。

`CLASSIFICATION_PLANE_M3_TRIPLE_VALID_ORE_MODULES_20261108.json` の
whole family `E1DX_1` は generic には不整合だが、
次の二 fixed wall は consistent であり、proper actual D factor を証明していない。

| τ | 完全 affine W-kernel | report の projective image 上界 | proper D factor |
|---:|---:|---:|---|
| 4/21 | 6 | 23 | 未証明 |
| 1/12 | 6 | 23 | 未証明 |

report は両方で `G_empty=false,G_closed=true` と明示している。
これらが実際に rank四既約 quintic Poisson point を持つ、又は持たないとはまだ主張しない。

`PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_G_BOUND.md` §2 の
「全 consistent module は proper factor、残 module は不整合、六源の全既約源は空」
という強い文言を、次に置き換えて読む。

**十八 module の全 complex parameter support は、空、actual proper factor、又は
次元二十四以下の全像上界により G 非優越である。
上記二 wall は次元上界だけの source であり、全点 E認識と空性は未解決。**

元 first-symbol 上界を安全に B≤9 としても、全 Q-source dimension 四、
完全 W-kernel 六、plane 選択六／projective scale 一により

`B+4+6+6−1≤24`。

τ fixed のため B≤8 を使えば report の二十三が得られるが、
今回の G 結論には安全な二十四で十分。
全 cohom／central coefficients は actual normalizer 逆で戻して数えており、
canonical module の W六だけを global source dimension と扱っていない。

corrected-V2 九十源の全 G 閉鎖、及び新 corrected fixed source の exact z-factor は保持する。
これら二 wall を含む pointwise E 被覆又は旧六源全点空性は受理対象にしない。
