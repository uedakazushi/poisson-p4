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
全74成分を個別同定し、全成分下界74、rank4下界68を記録し、上界・一般五点・全分類を未証明と表示する。
各checkpointは完了宣言ではない。mainへの保存後も、許可の再確認なく次の有効な数学へ進む。

負の法線重み三源と二超平面・三次型の四22成分も独立検証済みで、下界74/68へ更新した。
相異なる三H+既約Qの全配置とsmoothQ+既約C3全型は被覆済み。全nilrank2五型のprimitive一次像は15以下。
全rank1および非零定数modular型は非優越と証明済み。G-2Hの残余はA=c=0で共通因子なしの最高部。Eは未被覆nil源、1+1+3の純nil、特異Q+C3のrank4三源とrank3 pure nil十二型、他次数分割と既約非正規を処理する。
