# 五つの実 core 三次初期式: adjoint と五 IFT 補正

2026-11-04。一回の145係数 ring 計算。source 三次・四次解は新たに求めない。
前提 `CORE_CUBIC_ADJOINT_PREFLIGHT` の全2450恒等式と保存二次 source を使う。

base は四変数四次以下70単項式。五つの独立 dual 方向はそれぞれ core 二次以下15単項式。
元 scalar DAG の17329祖先を QQ で評価し、全98 unit 方程式の base 二次以下零と
dual core 一次以下零を確認する。それより上の残差を零と置かない。

選択五関数の三次出力は定数 adjoint lambda0 による source 三次補正を引く。
L1 で割って G選択関数の二次を得て、実対角行列 D5 により U2=-D5^-1*G選択,2 を得る。
残り五関数 rows63,64,68,73,74 の raw 原点勾配は全98方向で零。
従って真の四次出力は

    Fi4=Fi4(Y2)-lambda_i,1(C)*E98,3(Y2).

未知 source3 は adjoint により消去され、source4 は原点勾配零により現れない。
保存全九変数 source2 を五 dual 方向へ代入して mixed derivative を得た後、

    Ki3=(Fi4|U0+sum_j (partial_Uj Fi3)|U0*U2_j)/L1

を正確に実行する。Gi,2|core=0 により L2 は不要。
五つの U2 と五つの Ki3 の全 L1 除法を、それぞれ70係数で検算する。
各実関数・補正式・元残差を別添に固定する。

## 四次曲線上の厳密な非零 gate

既受理の nodal quartic 正規化 Phi:P1→Proj(k[core]/(f,g)) の四 binary quartics へ
Ki3 を代入し、binary degree12 section を全係数で求める。
その一つが非零なら、Ki3 は (f,g) に入らない。
実二関数曲面 S2 は初期環が k[core]/(f,g) である整域なので、対応する実 Ki も S2 内で非零。
従って全 core は高々一次元となる。
全 section が零の場合にはこの上界を主張せず、それ以上の次数へ自動的に進まない。

## 実 R による別個の下界 gate

既に source-bound された同じ R0 の実 R DVR 曲線は原点を通る。
E4 は D0 上有限で、二つの正規化の integral surjection の base change から
Spec(E4)→Spec(D0) は全射である。従って R の generic prime 上に E4 の prime がある。
その R0 prime は Nnorm 投影の norm 方程式を外れることが、実 -81600*T^4 係数で証明済み。
Npair hypersurface S/(L) はその Nnorm 自身なので、上の E4 prime では L は単位。
よってその prime は S/J、すなわち core 商上にあり、core の次元は少なくとも1。
この下界は形式的な tangent-cone 次数からではなく、同じ実 source と有限射から得る。
非零 binary section がある場合にだけ、上下を合わせて core 次元を正確に1とする。

全五 cubic と binary section の実結果は添付 report/evidence に保存する。
一次元と分かっても CM・H0・socle・E4 または D0 の深さは結論しない。
四 graph の ideal 冗長性、完全 syzygy、node gluing も別の未決問題である。
新 source jet 求解・GB・追加点・自動次数引き上げは行わない。

再現入口:

    .venv/bin/python code/ambient_B18_core_cubic_initial_20261104_v63.py --output work/B18-core-cubic-replay
