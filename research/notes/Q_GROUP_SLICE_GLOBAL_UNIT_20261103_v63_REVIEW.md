# 全109切断での17群方向単元：独立監査

2026年11月3日。frozen producerを再実行し、秒数以外一致した。
自前の既存sparse150 Lie engineから全289作用素と全109変数線形形式を再構成し、全係数が保存値と一致した。
各項の重み、108変数全ての厳密正重み、行和16・Lie和−16を確認した。
従ってdeterminantの重みは0、a以外の正変数に依存できず、QQ[a]で再計算した
163840*a³*(9a+4)/531441が **全切断上** の恒等式となる。

親D(a)開で単元なので、identityに接する17次元group chart×109切断から126canonical affine空間への写像は
全ての該当点でétaleとなる。全Poisson idealのGL同変性により、schemeへの制限も局所積である。
既存74変数の厳密消去を合わせ、canonical affine完成環は各R点の完成環に17 smooth variablesを加えたもの。

非zero tensorのスカラーGm商は滑らかなので、Rregular9が証明された点ではcanonical affine26、projective25が正則。
スカラー方向が選択17群方向だけに入ると仮定しておらず、任意の非Qa点でprojective完成環が
文字通りRhat[[16]]になるという追加主張はしない。全域一意なorbit sliceも主張しない。
code/q_group_slice_global_review_20261103_v63.py と同名 REVIEW report が独立証拠を保存する。
