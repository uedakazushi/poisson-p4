# Q_bの実共通29 base：étaleな事実と、full Dへの不足syzygy

2026-11-03。これは有効なbase projectionの事実を保存し、**それだけではfull DがfJ_Zとは証明できない**ことを明記する。
凍結済みCI4やquadratic-zero packetを変更しない。

## 同一の実29 target関数

B18 quadratic-image packetが選んだ共通28個のphysical target rowsとrow61を、そのままQ_bの三つの実sourceへ制限した。
同じtarget関数であり、sourceごとに別のprojectionを選んでいない。
全150 source微分から得た29次detはQQ(b)で

    CB1: −1728 b⁸(b²−5)/125,
    CB2: −442368 b⁹/25,
    CA: 8064 b/25。

b=1では全て単元なので、三つのcompleted smooth source/branch環は**同じ実target base** R=C[[29]]に同型となる。
b²=5の幾何的例外をこのprojectionだけで分類しない。

CBの二branchの全meetingをC2=R/Iと書くと、actual CB completed ringはcanonicalに

    A_CB=R×_(R/I)R = R⊕tI,  t=(1,−1), t²=1

である。実共通base関数はdiagonal Rを与える。
残る全physical121 target座標を各branchへ制限し、平均を引くとt b_jとなる。
これらは全ambient座標なのでb_jはIを生成するが、**多くは冗長**である。

## exact colon criterionと、無効な短縮の訂正

CAもR-chartなので、その同じodd target座標をw_j(R)と書ける。
CB2とCAのclean intersection N=(f)から

    w_j+b_j=f c_j

で、c_jの一つが単元である。全odd-generator vectorへ可逆R-row operationを施してc=e0とする。
この操作は冗長座標も保持しなければならない。
J_B=(b_j:j>0)とする。

R⊕tIの全presentationは、b-vectorの**全R-linear syzygies**とY_iY_j−b_i b_jで与えられる。
product relationsのCA制限はfJ_Bとf(f−2b0)である。
一方、syzygy a0 b0+Σ_(j>0) a_j b_j=0の制限はf a0。
可能なa0のidealは正確に(J_B:b0)だから、actual D idealの全recipeは

    Q_D=f*( J_B, f−2b0, (J_B:b0) )。

CB1とCAのactual intersectionはJ_Z=(J_B,f−2b0)なので、product relationsだけからQ_D=fJ_Zとするのは不十分。
common baseがétaleでも、冗長ambient座標のsyzygyがf²rhoを追加し得る。
conductor Iがregular sequenceで最小9generatorを持つことだけでは、**全実ambient odd vectorのc方向**を最小9へ落としてよいとは限らない。

相対CI4と実isotropic quotientにより、Z4ではf²=rho²=0、ann(f)=(f)である。
従ってcolonのZ4像は(f)に入り、先行actual quadratic-zero gateはそのlinear f classを除く。
残る可能性はsocle f*rhoの一方向、すなわちactual Dへのcubic f²rhoの追加である。
この全121-coordinate bitを実third-order source inversionで確認することが次の課題であり、ここでは未解決のまま保存する。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CB_self_pair_common_base_etale_20261103_v63.py --output <new-output-directory>

全projection行番号と三つの厳密detを保存する。GBなし、120秒・2.2GB。
