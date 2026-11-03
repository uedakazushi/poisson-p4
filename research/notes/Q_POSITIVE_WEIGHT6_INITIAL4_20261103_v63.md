# 正 weight6 軌道の全 ordinary initial 三次・四次と局所次元上界25

2026-11-03。実八変数形式模型の**真の ordinary initial ideal の全次数4まで**を求めた。
二次5本に対し、新しい三次最小類が8本、四次最小類が8本ある。
実横断局所環の Hilbert function は次数4まで厳密に

    (1,8,31,72,126)

である。これら21本の斉次式で生成する有限 polynomial ideal の次元は3、次数は6であり、
従って実形式横断環の次元は3以下、全 projective Poisson 局所次元は25以下となる。
非零正軸の全点は一つの既証軌道なので、上界は全非零軸に適用する。
実次元3、全 initial ideal の生成、実枝の存在・本数・尽くしは主張しない。
特に次元24以下という次の目標はまだ証明していない。

## 全350式の三次・四次解

凍結済み `Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63` の23 orbit cut、
95元式の単元線形 block B、八方向 T をそのまま使用する。
canonical raw tensor は affine 線形切断のため、形式解は

    P(z)=P0+R1(z)+R2(z)+R3(z)+R4(z)+...

と書ける。R1=Tz、Rdは dependent95方向だけを使う斉次次数dの補正である。
strict Jacobi map は物理的 tensor に関して二次なので、degree3の未補正項は
Bilinear(R1,R2)+Bilinear(R2,R1)、degree4は
Bilinear(R1,R3)+Bilinear(R2,R2)+Bilinear(R3,R1) である。
各次数の選択95行へ −B^(-1) を掛けて Rd を一意に定め、
**全350式**へ戻して残余 Q2,Q3,Q4 を記録した。選択95行の三次・四次残余は全係数で零となる。

全八変数の QQ 係数を保持した有限演算であり、標本や modp の階数から三次・四次式を推測していない。
追加検算として、三つの dense QQ 八方向で長さ5の scalar-series tensor を別に組み立て、
物理的 Jacobi の全350式を univariate convolution で直接計算し、次数0–4の全残余と照合した。
この標本検算は全係数証明の代用ではなく、独立な indexing・sign・次数の検査である。

## degree3 cancellation の完全性

Q2 の独立な元式五行を A=(A1,...,A5) とし、全残余350式からその定数線形結合を引く。
選択五行を除いた345本 G はすべて次数3以上から始まる。
これは形式級数生成 ideal を保存する可逆な定数 row 操作である。

全 degree3 initial 空間は

    S1*A2 + span(G3)

である。5*8+345=385行、120列の完全な QQ 行列を保存した。
rankは48、したがって h3=120−48=72である。
40行の線形倍を除いて新しい8本の cubic class を選んだ。

## degree4 cancellation の完全性

degree4の新しい initial を生むすべての低次 cancellation は、385行の上記写像の kernel に入る。
その全 QQ kernel は次元337で、全 kernel 行を元行列へ掛けた零性を検証した。
各 kernel の同じ係数を

    (S1*A3, G4)

へ掛けると degree4 carry を得る。全 degree4 initial 空間は正確に

    S2*A2 + S1*G3 + carry(kernel)

である。最初の二項は下次数の initial の倍で、最後が cancellation の全補正である。
rankは204、従って h4=330−204=126。下次数の倍を除いた新しい quartic class は8本となる。
各元式の最初の非零次数だけを集めた計算ではなく、全337 cancellation を含む。

完全性は次数を数えるだけでも確認できる。A2の五行は定数線形独立なので、
degree4で始まる任意の ideal 元の A への定数係数は零である。
次の次数で消える条件は Aへの線形係数とGへの定数係数が上記 kernel に入ることであり、
残る Aへの二次係数、Gへの一次係数は最初の二項を与える。ほかの degree4 source は存在しない。

## 次元上界と保存

小さな QQ 八変数 ideal の正確な GB により、二次と三次だけでは次元5、次数1、
全四次まで加えると次元3、次数6となる。実 full ordinary initial ideal はこの有限 ideal を含むので、
実形式環の次元は3以下である。23 affine smooth因子、または22 projective smooth因子を戻し、
全 projective 局所次元25以下を得る。逆包含や次元の等号は使用しない。

producer は `code/q_positive_weight6_initial4_20261103_v63.py`。
`reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63/` に、全350三次・四次、95形式補正、
385行の全写像と337次元kernel、全 carry、21本の generators、QQ GB入力とログを保存する。
最初の全係数計算は約28秒で終わり、大規模 GB や新しい五次 jet は実行していない。
