# 正weight6軸：全350の五次係数と全initial5の独立検証

2026-11-03。凍結済みQ_POSITIVE_WEIGHT6_INITIAL5 packetを検証する。
前回reviewのcanonical126・orbit23・implicit95・形式8表示を再使用し、
物理的QQ全係数と五次の全低次cancellationを独立に確認する。
115秒のproducer全計算は再実行しない。

## 1. 標本でない全QQ物理係数の検査

reviewerの五物理変数から作った二次積表と別C++/FLINT engineを拡張し、
P0+R1+…+R5を八変数の斉次次数に分けて厳密に畳み込む。
全350式×次数0–5の2100個の多項式、divergence25×6次数を比較する。
全八変数QQ係数を保持し、方向標本やmodpの検査ではない。
選択95元式の五次残余が全係数で零であることも確認する。

## 2. polynomial row eliminationと全五次source

最小二次A5、最小三次C8を使うdegree3の48行は独立である。
元G345の三次部分をその48行で全て表示し、同じ係数を真の形式式へ掛けて引く。
こうしてorder≥4のH337を得る操作は可逆な多項式row操作で、形式idealを保存する。
全四次部分と五次部分の同じrow表示をQQで再計算する。

最小四次D8を選ぶと、degree4 sourceは

    M4=(S2*A2, S1*C3, D4), 252行、rank204。

保存48行kernelについて零積とrank48を確認し、252−204=48から全kernelと分かる。
同じkernelを(S2*A3,S1*C4,D5)へ掛ける全48個のfifth carryを再計算する。

H337の四次を204独立行で表示し、同じ係数で元形式式を引けば、order≥5のE329を得る。
その337×204表示、全四次消去、五次残余329本をそれぞれ検査する。
従って形式idealはA5,C8,D8,E329で生成される。

degree5 sourceは正確に

    S3*A2 + S2*C3 + S1*D4 + carry(kernel M4) + E5

であり、952+48+329=1329行、792列となる。
degree2の定数A2、degree3の(S1A2,C3)にはkernelがないため、
それ以前の未計算cancellationは残らない。
degree4の全kernelとorder5 generatorを含むので、五次のsourceはこれで尽くされる。

## 3. exactQQ rankと二つのlinear support

有限体で選んだ新五次15行が、保存された実1329行の指定行そのものであることを確認する。
独自に作った小さいQQ八変数GB入力で、残り全377行(E329+carry48)が零へreduceすることを検査する。
下次数の倍は生成idealに元から入り、逆に選択15式と下次数は実initialなので、
この包含はdegree5空間の等号を与える。
したがってMacaulay2のhilbertFunction(5)=171から、真のQQ degree5 rankは792−171=621と確定する。
modp rankをQQ上界へ外挿していない。
実Hilbert functionは次数5まで厳密に(1,8,31,72,126,171)である。

二次radicalの二つの5平面へ制限し、そのQQ radical idealを独立に比較する。
最初の平面では(u6,u4,u2*u3)、第二では(v6,v2*v3,v2*v5)となる。
全五次coneの被約支持は、包含を整理すると

    L3: free(z0,z3,z5),
    P2: free(z0,z2)

の和集合、すなわちradical=(z1,z4,z6,z7,z2*z3,z2*z5)である。
有限idealの次元3・次数5を確認するが、実初期ideal全体や実局所環の次元・次数と同一視しない。

## 4. 結論の範囲

全projective Poisson局所次元の上界25は維持される。
五次まで残る3平面が真の接錐に全て含まれること、実3次元枝が存在すること、
全枝・全初期idealの生成を示したわけではない。
次の数学的課題は、この3平面の実際の持上げ、または高次obstructionによる切断である。
正軸のprojective安定化Lie代数はGa²であり、点を固定しないeta移動の対角作用を
正のgradingとして使うことはできない。

review code: `code/q_positive_weight6_initial5_review_20261103_v63.py` と同名`.cpp`。
receipt: `reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW/`。
