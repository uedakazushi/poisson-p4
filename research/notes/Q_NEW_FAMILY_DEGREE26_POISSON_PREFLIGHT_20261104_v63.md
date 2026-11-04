# 26次の普遍基底をPoisson括弧で構成する案：事前検査

2026-11-04。全Inv₂₄の排除が独立受理され、NEW25の定義次数は26以上となった。本packetは**同じ固定形式ambient平面を用いる、全Inv₂₆基底構成の有限な方法**を検討するだけで、26次のrank計算・source展開・新しい点評価を行わない。

既存のexact Weyl character計算から

\[
 \dim\operatorname{Inv}_{14}=3,\quad
 \dim\operatorname{Inv}_{18}=9,\quad
 \dim\operatorname{Inv}_{26}=45
\]

で、25,27次は零。rootが保存した全Weyl8項の記録を26次まで同梱した。traceが全不変式環を生成するとは仮定しない。

## H₁₆の微分を避ける候補

受理済みの五つの4×4 C/D trace4

\[
 (0000,0001,0011,0101,0111)
\]

は全Inv₈の基底であり、canonical pᵢの基底との変換も既に証明されている。その全10個のPoisson括弧はInv₁₄を張る。旧Inv₁₄の3次元基底も全Inv₈の括弧から得ているので、basis changeによってimageは変わらない。同じ形式平面の係数でもrank3となることを後続gateに要求する。

12次のC/D候補はword000101、000111、011011の三つ。最後は既存N₁₂に対し

\[
 N_{12}=T_{011011}-F_4T_{0011}/20+F_4p_4/5
\]

を満たす。一方、任意のSp₄不変式fについて

\[
 \{F_4,f\}=0.
\]

実際、moment mapの定義からF₄のHamiltonian vectorは、moment係数を持つ10個の物理的な無限小Sp₄作用の線形結合となるため、不変式を殺す。従ってN₁₂とT₀₁₁₀₁₁を交換したときの括弧の差はF₄Inv₁₄に入る。

よって次の18候補を同じ形式平面で検査する案は、N₁₂を使ったrootの案と同じspanである。

* F₄Inv₁₄の3本。
* 五つのInv₈と上の三つのtrace6のPoisson括弧15本。

最初の3本を含むrank9が得られた場合だけ、その補空間の6本をprimitive18として選ぶ。これにより

\[
 3\operatorname{Inv}_{14}\cdot9\operatorname{Inv}_{12}\quad(27本),
 \qquad 6\operatorname{primitive}_{18}\cdot5\operatorname{Inv}_{8}
 \quad(30本)
\]

という57本の26次候補を作り、rank45を要求する。rank不足を候補全体の普遍的上界と扱わず、また自動的にH₁₆のgradientへ拡大しない。

## ambient16変数の微分と形式計算費用

**C,Dの全136 quadratic係数から周囲の16個のβ座標で微分した後に、固定二平面またはsourceを代入する。** s,tだけの微分ではΩ⁻¹によるPoisson括弧を得られない。source上で零になっているβ座標の微分方向も捨てない。

trace wordの微分にはcyclic product ruleを使う。全必要prefixは長さ1,…,5で2,4,8,13,14個、4×4 scalar polynomial product上界2496回。8個のtrace invariantの全16 gradient contractionは3840回以下。Ω⁻¹の非零成分は24個で、10+15個のPoisson括弧は600回のscalar productで計算できる。符号は元の規約∇fᵀΩ⁻¹∇gを保持する。

形式二平面ではgradientの次数は7または11、係数数はそれぞれ36、78。全gradientでも6624 residue程度、57本の26次係数行列は378×57。rank gateは3、9、45。全入力、全semantic word、係数配列と非零minorを保存するnative一回60秒/2GB/certificate32MiBという案は小さな有限計算となる。まだ実行していない。

## 将来source計算の分母・奇偶性と、上限保証の不足

sourceでは既存のT²=Bによる厳密なtorus前処理を再利用できる。4C,4Dのambient gradientは整数、5Ω⁻¹も整数なので、Poisson括弧の整数分母はcontent除去前に

\[
 \operatorname{Inv}_{14}:\quad5\cdot4^8\cdot1080^{14}B^7H^{14},
\]
\[
 \operatorname{Inv}_{18}:\quad5\cdot4^{10}\cdot1080^{18}B^9H^{18}
\]

と追跡できる。26次の全積はB¹³H²⁶を共通分母に取れる。

cのhomogeneity次数は14,18,26次で7,9,13となる。既存の有限τは(u,X,c)を同時に反転するため、c=1後のsource式は**(u,X)に関し奇**である。従ってu=aXを先に代入すると共通X因子を持ち、それを除けばY=X²の多項式になる。u,Xの多項式段階でXが因子だとは主張しない。dense open X≠0上の有理変換なので恒等式検査は保たれる。T次数は全て偶数となり、B=T²へ戻せる。

既存整数source numeratorのℓ¹ normと全16方向のderivative係数から、**新しいgradient多項式を展開せず**次の厳密なmajorantを得た。

|次数|X除去後support上界|整数係数bit上界|出力一つのraw bytes上界|
|---|---:|---:|---:|
|14|48420|556|10071360|
|18|124730|712|31806150|
|26の列|508760|未展開|未保証|

必要な新source mapsはInv₁₄三本とprimitive18六本。最悪raw上界は約221MBで32MiBに収まらず、18次supportの理論上界も100k項を超える。したがって**source計算の完走はこのpreflightでは保証しない**。まず普遍基底gateと選択wordを確定し、別の明示的budgetで60秒/2GB/出力100k項/中間150k項/全raw32MiBなどの停止条件を再検討する必要がある。実測が上界より小さかった前回を今回の保証に使わない。

保存：`code/q_new_family_degree26_poisson_preflight_20261104_v63.py`、
`reports/Q_NEW_FAMILY_DEGREE26_POISSON_PREFLIGHT_20261104_v63/`。
全character26次までのWeyl記録、word・prefix定義、16方向のmajorant、分母・parity、input hashを含む。metadata算術は約0.22秒。普遍26次基底、source関係式、次数26の方程式・排除は本packetでは主張しない。
