# 次の課題 — 全被約成分の尽くし

2026-11-07 / v64-classification-first。classification/PLAN.jsonが現行キュー。

```sh
.venv/bin/python scripts/classification.py check
.venv/bin/python scripts/classification.py next
```

G、E_multiple_hyperplane、E_reduced_nonnormalを並行する。
Gの入口はclassification/DOMINANT_COMPONENT_GEOMETRY.mdとREMAINING_DEGENERACY_DIVISOR_SUPPORT.md。
未知beta優越成分は射影次元25〜27。正規な実退化因子の全点は既存91頁原稿の定理でT39へ分類済み。
これを見落としてG-normalを未証明とした旧草稿のキューは撤回した。
局所Lie補題は補足として保存し、閉じた正規分岐の再分類は行わない。

全134式のnormal rank4枝とF5、有限分離忠実chart、直接splitterは利用可能だが、
それらだけからB=0とはしない。独立reviewの五点+孤立第六点反例が限界を示す。
五回未完了の直接GB、低次trace探索、旧compact方向graphの再包装を反復しない。

E_quadricとE_quadric3は完了。実退化因子H+2Q、Q既約(rank>=3)の全支持を覆った。
新22次元混合型と新24次元Tw0型を全周囲の真正初期障害と実源下界から証明し、独立再構成も完了。
正規因子と二重既約二次因子を持つ全支持を再分類せず、その補集合へ進む。
可約Qは二重超平面を持つ側へ残す。

E_multiple_hyperplaneは全affine二次P=J+L+Q、全Jacobi、Pf(Q)=0を保持する。
modular場を使って全固有値/Jordan分岐を有限に扱い、既知Q45結果とSNC実変形を適用する。
E_reduced_nonnormalは被約非正規五次因子の全支持を扱い、既知MPS成分を上界と仮定しない。
新しい非自明な計算は該当active taskへの有限job契約を作ってから行う。

E_zeroとE_poly、A_rank2は完了。従前58成分のA監査と新22/24次元成分の独立証明も完了。
旧証明書を入力にする範囲と今回の独立再計算を区別し、新候補だけを必要に応じて追加監査する。
全86成分を個別同定し、全成分下界86、rank4下界80を記録する。上界・一般五点・全分類は未証明。
各checkpointは完了宣言ではない。mainへの保存後も、許可の再確認なく次の有効な数学へ進む。

既約Q＋既約C3の全支持は全cofactor型を尽くし、MPS40と六つの新23成分へ被覆済み。
平面二直線＋既約三次は十三全配置へ尽くし、六23成分を証明、残七配置の帰属を保持する。
nil31 long-top六源とnil22 top八源も全被覆。short r=2/3、image内nil源等は未被覆である。
全rank1、非零定数modular、全nilrank2型はprimitive一次像を優越しない。
GはR13/R22/S23最高部の全閉包を非優越へ送り、logarithmic/exceptional最高部と全境界を処理する。
Eは上の残支持、H＋既約quartic、他次数分割と既約非正規quinticを継続する。
