# 真正四次元 target の E7 と kernel：独立値検算

主報告 `c93b6250...` に対する独立報告 `f35d75c4...` を固定した。
旧 C2 の24葉を出典から読み直し、四演算子の全60項を逆順に加算した。
全 termwise product と四値

    alpha_A=-85/6, alpha_B=8315/432,
    beta_A=275/3, beta_B=-359975/1728

が一致した。新しい source、physical、gauge、IFT、jet は使っていない。

実 E7 と7×4 kernel行列を独立構成し、全28積、rank三の単位 minor、
kernel の自由 I4を確認した。E²の像は非零 beta vector による一列だけで、
E³=0、rank列は7、3、1、0である。
従って独立 SOURCE 監査済みの同型を通じ、真の C[ell] 四次元の座標が得られる。

これは選んだ primitive D-free F に関する D 商の座標である。
F は R0 安定でないため、この E7 から素朴に r/z 作用を降ろさない。
R0 transport、全形式 torsion lift、connecting 値、元 H0 や nil 三乗は別の課題である。
再現は `code/ambient_B18_C_genuine_T7_values_review_20261105_v63.py`。
