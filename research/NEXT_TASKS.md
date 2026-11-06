# 次の課題 — 全被約成分の尽くし

2026-11-07 / v64-classification-first。classification/PLAN.jsonが現行キュー。

```sh
.venv/bin/python scripts/classification.py check
.venv/bin/python scripts/classification.py next
```

主課題Gは未解決。まずclassification/NORMAL_DIRECTION_INCIDENCE.mdと独立reviewを読む。
全K16で交換子[T0,T1]が可逆であることを非零minor21866で示し、
方向incidenceのrank3枝を除外した。残るのはrank[A b]=rank[A b d]=4上の全F5である。
明示的な多項式左逆から既知三平面と未知残余のsplitterも直接得られるが、残余Bの空性は未証明。
この残るrank4方向集合を全て処理するか、元134式を満たす未知解・族を構成する。

classification/FINITE_SEPARABLE_NORMAL_CHART.mdにより、K5へ有限48/144被覆を施せば、
42独立係数を持つ一つのellが未知残余を忠実に検出する。追加の超越変数は不要。
全134式・259法線式への引戻しとGalois指標による根基の復元は証明済みだが、unitは未取得。
一つの数値根だけを選ぶ特殊化へ置き換えない。分岐補集合と全底分母をEへ残す。
旧全K16固定etaのK5への単純コピーは、この別形式の定理では正当化されない。

直接GBは過去五回未完了。normal-only展開も今回四次式前に331,008項と判明した。
同じ点・素数・予算の直接GB反復を行わず、未展開の行列構造または支持幾何を使う。
計算前には新しい命題と分類への矢印、全補集合、有限契約を明記する。

E_zeroとE_polyは処理済みの明示部分である。
beta=0は旧Q45_22を入力に既知52閉包で覆われる。
SNC最高部を持つ二次多項式代表では全Poisson低次項を覆えるが、
Qの閉包上だけにある垂直成分、他の最高部、三次尾部を持つ全構造は未処理。
Gの底分母は既知49proper像で消え、G28を因子に持つ必要がある。十分条件とはしない。

A_rank2は完了。RANK_TWO_SUPPORT.mdに全安定・非安定・因子付き分岐、
CM/LPT一次本文の適用、全Aff32変形、独立査読がまとまっている。
六全空間rank2成分の尽くしは、rank4の全列挙なしで成立する。
Aの残りは追加11rank4族とNEW25の全周囲成分性・相違の証拠を主原稿へ接続すること。

主原稿papers/poisson_p4_classification_ja.texはこれらの証明を統合している。
全成分下界58、rank4下界52は維持し、上界・一般五点・全分類を未証明と表示する。
checkpointは完了宣言ではなく、常駐研究も開始しない。
