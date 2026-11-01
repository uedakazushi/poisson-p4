# free-M0 八次有限上界の証明・source監査

2026年11月1日。[主証明](S61_MINIMAL_FREE_FLAG_EIGHTH_UPPER_BOUND_20261101.md)
とraw8 producer、切断standard-basis nativeを読み、関数体上の
`length(S/m9)≤6745`、`HF7+HF8≤4177`を受理した。
独立のraw8・standard basis再実装や主producer再実行は行っていない。
入力hashと整数範囲・次数の算術だけを別checkerで再確認した。

nativeの各追加行は、元の切断式および既に構成した行への単項式乗算、引算、
非零有限体係数による除算から得られる。従って各行は元式の有限Macaulay行空間に属する。
これはS-pairの全処理やproduct/chain criterionの完全性を仮定しなくても成り立つ。
各先頭単項式の倍のうち次数8以下の異なる68837個について行を一つずつ選べば、
先頭が互いに異なるので三角性から独立である。
この監査は68837をrank下界として使い、有限体の全rankとの等号を独立に証明するものではない。

元135式の定数項とJacobianの照合、陰関数再帰のp単元分母、
全項を保った八次までの計算はsourceで確認した。
先頭式だけを元の正確な方程式へ追加する操作は行っていない。
整数内積とbincountの和はそれぞれint64と厳密な倍精度整数範囲内である。

有限行列の68837次非零minorはQQ(a0,…,a4)でも非零であり、
75582−68837=6745が一般m9商長の上界となる。
真六次までの一般総長2568は先行の普遍定理で確定しているので、差4177を取れる。
一般HF7が確定していない以上、指定点の差2410を一般HF8の上界にすることはできない。
主証明がこの限定を維持していることを確認した。

全初期ideal、tの冪零性、飽和除去についての結論はない。

```
.venv/bin/python code/s61_minimal_free_flag_eighth_upper_bound_review_20261101.py --output /tmp/free-flag-eighth-upper-review.json
```
