# 正 weight6 軸の既知成分排除の独立監査

2026-11-03。`Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63` の結論を承認する。
a=8/9 の実曲線の全 eta≠0 点は actual T01、CA、CB、および共通被約 N の外にある。
他の既知成分への所属、新成分、局所枝本数、全模型の非平坦性は結論しない。

物理的 tensor の六つの成分を別に組み立て、保存された全150係数と一致させた。
divergence25係数と strict affine Jacobi 自体を直接微分して零と確認したため、射影350係数も零である。
実五次式は −16*x0^5/81 である。固定 Darboux 変換後の全61座標、perspective 零、
全150復元・元の134式を検証した。

閉二次 jet の規約は物理的 jet から独立に確認した。Q を二次 bivector とし、
二形式 −JQJ の Euler contraction を外微分して4で割り、再び J で上げる。
この直接計算は E*Pi*Q と一致する。したがって保存式の入力は Pi であり、
分母払いの producer 内部の4Piではない。一次 jet の primitive 補正 alpha は零である。
403項の保存 highest relation を係数辞書として読み、別の積・和演算で評価した。
生き残る四項の和は厳密に −1163264*eta²/135 となる。

元 producer を新規 work ディレクトリで再実行し、全274式の QQ 根作用安定性814像と
六変数 coframe 全係数恒等式を含む report・証明書が凍結物と一致した。
全 T01 に及ぶ論理も確認した。六変数 R0 は交代形式 RJR^t の稠密開の代表であり、
同じ形式の coframe は R0*Sp4 で尽くされる。四つの根作用で安定な274式の空間は
連結 Sp4 に安定なので、稠密 coframe での零性は actual T01 の被約閉包全体へ延びる。
これは discriminant で割った識別式ではない。

coframe jet の分母払いの係数も式から点検した。
F=R^(-1)(I−Ry*1^t/5)D R、D=diag(1+Ry−sum(Ry)/5) の一次係数を A_i とすると、
producer の行列は5*det(R)*A_iである。したがって線形 bivector は10*det(R)倍となる。
二次の FJF^t を展開し、4Piを掛けた producer の値は
(10*det(R))²*Pi(Q) である。重み4の恒等式の零性はこの非零共通 scaling で変わらない。
単一の単位 coframe の標本較正だけを普遍的較正の証明にしていない。

CA/CB の排除では、物理的 primitive beta=−JLJ から Kbeta を直接再計算した。
三次 minor は −384*eta で、eta≠0 上の全射影 kernel は唯一 e_w である。
固定 Sp4 変換はこの旗を reference e_y に送り、実61座標への有限作用後、
両方の完全な閉 beta-image 制約に非零定数3/2が現れる。
既承認の proper closed flag-source は S61 の全被約 CA/CB を覆い、全候補旗はこの kernel 内にあるので、
未検査 chart の旗や別の旗に逃れる可能性は残らない。

独立コードは `code/q_positive_weight6_exclusions_review_20261103_v63.py`、
receipt は `reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63_REVIEW.json`。
新しい排除式は全非零軸に適用するが、ほかの成分への所属は引き続き未解決である。
