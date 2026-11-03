# 真の四次初期式の全377 cancellation 被覆の独立監査

2026-11-03。完了済み202秒の full DAG 四次 jet・独立 physical350再検証の source と report を
監査し hash 固定した。理由なく同じ大計算は繰り返していない。
保存 raw350 四次 tail を用いた bounded mixed producer を新規再実行し、
展開済み証明書、native整数入力とkernel出力は一致した。

元350式から57個の選択 original A式と293個の G式への定数行操作を独立に追跡した。
各 G は自身の元行を係数1で含み、他は A57 の組合せなので全取替えは可逆である。
全293 cubic tail は既監査の三次資料と一致した。
219個の recipe の三次恒等零と actual350 の定数＋一次 cofactor を独立に確認した。
さらに weight block ごとの mod1000003 rank 下界219と、既存の exact QQ image rank1911から、
この219方向は1995+135三次行の kernel 全体である。
零 cubic の158座標方向を加えた377方向が full kernel を漏れなく張る。
101線形 syzygy と55定数 cubic syzygy を219方向へ重ねて数えてはいない。

両15平面では、保存済み reduced quartic carry の値に頼らず、
元の raw350 quartic D_i を QQ[15] に直接制限し、actual350 cofactor の定数項×D_i と
一次項×C_i を掛け直した。全219方向と零 cubic158方向の carry が全て零となった。
全350Q と293Gの cubic 制限も直接零を確認した。

任意の ord≥4 の I の元を A57,G293 で書くと、独立な Q57 により A の定数係数は零である。
三次 cancellation は正確に377次元 kernel、四次項の残りは S2Q と S1C_G だけとなる。
従って (in I)_4 は S2Q、S1C_G、kernel 全 carry の spanであり、今回の平面上の消失は
真の四次 component 全体について成立する。
四次以下の全初期式の ideal は両15平面を含み、任意12本の同次線形切断の次元は3以上である。
full initial ideal・実模型の次元下界15や実15次元枝は結論していない。

監査結果: `reports/Q_ORDINARY_QUARTIC_INITIAL_20261103_v63_REVIEW.json`。
