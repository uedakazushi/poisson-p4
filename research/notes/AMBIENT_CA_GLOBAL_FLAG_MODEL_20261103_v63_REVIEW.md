# 全射影 CA の proper 旗模型と階数二境界：独立監査

2026年11月3日。新しい作業ディレクトリで frozen producer の全 QQ 証明書を再実行し、
実行秒数以外一致した。追加コード ambient_CA_global_flag_review_20261103_v63.py は
添字・単項式に対する独立の sparse150作用で、全20基底の canonical 化と320個の Lie 作用列、
二つの境界の全微分行列、小行列式、175×5 の定ベクトル因子行列を作り直す。
全20パラメータの生 Jacobi と実退化五次式も独立に展開した。

旗は L⊂V の直線、W⊂ann(L) の二平面、H⊂W の直線で、次元は4+4+1=9である。
同じく annihilator の旗 L⊂ann(W)⊂ann(H) から F(1,3,4;5) と見られる。
GL5 旗安定化群は連結 parabolic 次元16である。
その Lie 作用が保存する canonical20空間は群全体でも保存され、随伴 vector bundle が定まる。
固定 raw 空間の全0i成分は零なので、Euler-vertical な元なら hAi−XiA0=0 から A=cE、
従って raw tensor 自身が零となる。よって canonical 化は単射であり、projective bundle から
全 canonical Poisson 空間への写像に basepoint はない。
滑らかな projective28次元源の像は proper・閉・既約で、g≠0 の dense 開が既存 CA の
g=1 族そのものなので、像は全被約 CA と一致する。

実 rank≤2 条件は g*h³*f≡0 と同値であり、全 source の被約逆像は g=0 または f=0 である。
二つの線形部分空間について全16作用行列の対応 subblock が零であることを追加確認した。
従ってそれぞれ全旗上の projective subbundle を作り、源の次元は27と25である。
保存された二つの actual150 微分行列を独立に再構築し、非零 QQ minor のサイズ28と26を確認した。
affine image cone の次元下界から projective 像の次元は27と25で、源の上界に等しい。

条件「ある非零定ベクトル v に対し E∧v∧P=0」は正確に175×5行列の rank≤4であり、
代表を Euler-vertical に変えても不変、GL5 下でも不変な閉条件である。
g=0 の全 source は v=∂y を持つ。f=0 の保存点では独立に再構築した行列の rank が5である。
従って小さい25次元像は大きい27次元像に含まれない。逆包含は次元から不可能である。
二つの proper 既約像が全 boundary を覆うことと合わせ、被約 rank≤2 境界は正確にこの二成分である。

この結果は S61 に制限しない全射影 CA の境界支持の結果である。
源の smoothness から target の normality・R1を推論せず、境界における旗ファイバーや CB との
所属を結論しない。全 Poisson 空間の成分数へ新しい成分を加えるものでもない。
証明書は reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json に保存する。
