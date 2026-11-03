# CA/CB/NのR1と指定Qの非CM性の独立監査

2026-11-03。固定六平面のprimitive二形式を明示係数から別に構成し、
全旗方程式、11×3の縦接行列、その全最大小行列を保存証明物と比較した。
旧common_plane関数を読み込まず、全多項式が一致した。
独立Cloud M2は、最大小行列のideal Hについてt³∈Hと
H+(t)=(t,ΔA,ΔB,ΔC,ΔD,ΔF)、およびH⊂(t,Δ)を確認した。
Δ=4CF−D²なので、critical supportは正確に(t,Δ)である。

微分のcritical locusだけで全旗ファイバーが一点とは結論できない。
節点の正規化のような二つのunramified preimageの反例を排除する必要がある。
ここでは本文の全射影旗に対する直接論証を確認した。
t≠0ならKβvからv_x=0となり、chart v_y=1の実際の収縮式が
(Aq+Bw+2t)q/w、および(Aq+Bw+t)q²/w²を与える。
この四式はq=w=0を強制し、chart外は収縮のt w,−t qで排除される。
t=0,Δ≠0では、もしv_x≠0なら零のy成分から比例定数が零となり、
収縮のHessian block v_x[[2C,D],[D,2F]]が零となって可逆性に反する。
従ってこの場合も基準線だけである。接行列のrank3によりその一点は被約である。

この全ファイバーの議論があるため、源のcodimension2 critical集合のproper像の外では、
各成分解消のファイバーは空でない被約一点となる。
proper+quasi-finiteから有限であり、B/m_A B=CとNakayamaからB=Aとなる。
従ってその外で解消は同型で、滑らかな源の正則性をtargetへ移せる。
悪い像の次元は源のcritical集合の次元以下なので、被約CA/CB/Nは全てR1を満たす。

指定Q（a=1）では、既に独立に確認した有限正規化ファイバーが長さ2であるため非正規である。
もし対象局所環がCMならS2を満たし、R1と合わせて正規となって矛盾する。
従って三つの被約成分自体の局所環は非CMである。
それぞれの完成環がdomainで重複度2であることは、前の正規化監査から保たれる。
通常cuspと滑らかな空間の積はhypersurfaceでCMなので、その局所模型も排除される。

これはa=2/9の全35変数模型に対する別の非CM証明と、対象を混同しない。
全局所生成関係、conductor、全Jacobiの枝や埋込み構造は未決定である。
再現器は `code/ambient_component_r1_review_20261103_v63.py`、
証明物は `reports/AMBIENT_COMPONENT_R1_NONCM_20261103_v63_REVIEW/`。
