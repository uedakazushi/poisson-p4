# CA 原点 self-flag 三次 jet の独立監査

2026-11-03。`AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63` の全150係数と有限 witness を、**有限5×5行列の物理的 pushforward** から独立に検証した。作者の Lie / logarithm / jet helper は import しない。この監査は ordinary chart の有限3jetを受理するものであり、全次数の形式 branch や conductor を結論しない。

## 1. 全 source と旗方向

raw CA 式を canonical 化して固定旗20基底を自分で再構成し、受理済み critical 証明書の全150×20行列と照合する。実際の P0 も `∂y∧D0` から再構成する。既受理の parabolic orbit10方向を取り除いた first source の自由座標は `[1,2,3,9,10,14,15,17,18,19]` である。

relative flag の全9方向は `[E02,E03,E04,E10,E12,E13,E14,E32,E42]`。その全一次 quotient の階数は8で、核は

\[
M=E_{10}-2E_{03}
\]

の一方向である。ordinary separation では E10座標を parameter aそのものにできるので、残る8方向すべてを B,C3 に残し、

\[
g=1+aM+a^2B+a^3C_3
\]

とする。この opposite-unipotent chart の一般行列 N は N⁴=0、det(1+N)=1を満たすことも全9記号で確認する。

## 2. 別の coefficient 計算

独立コード `code/ambient_CA_origin_self_flag_3jet_review_20261103_v63.py` は

\[
g^{-1}=1-aM+a^2(M^2-B)
 +a^3(-C_3+MB+BM-M^3)\pmod {a^4}
\]

を使い、`g P(g⁻¹X) gᵀ` の各 tensor 単項式を直接変換する。二つの出力 wedge 因子と二つの入力単項式因子の積を a次数3まで厳密に計算する。`L_M²` と `L_(M²)` の間の同一視、log展開、作者 helper は用いない。

固定旗20空間の正確な左逆から作る quotient projectionを、全150行に適用する。二次にはP1/free10とB/flag8の全18変数が現れ、その150×18行列は保存値と一致する。選択15 minor の逆と全列 factorizationにより階数は厳密に15、consistentで、核は3次元となる。自由列は `[12,14,17]`、すなわち `B04,B13,B42` である。

三次には高次 first-source P2/free10 と flag C3/8が同じ18列の image で現れる。これは a次数を数えても直接の物理計算でも確かめられる。低次変数の全99個の係数列を直接作り、全150×99の raw quotient と、全18高次変数を消した cokernel を保存値と照合する。Bの二つの因子はa⁴以上なので、三次でBの基底方向だけを評価することは厳密な全列計算である。

## 3. 全 cokernel ideal

二次の全 affine 解を3変数vで代入した全150行は、独立な三行

\[
q_1=\frac{2v_2(2v_0+3v_1+5)}5,\quad
q_2=-\frac{4v_0v_1+4v_0+7v_1^2+18v_1+8}2,\quad
q_3=-\frac{v_2(v_0+4v_1+5)}5
\]

と同じ QQ定数 row spanを持つことを確認する。従って少数の選択 rowだけの検査ではなく、全150行の idealである。

`F=4v0v1+4v0+7v1²+18v1+8`、`H=4v0+7v1+11` とすると

\[
v_2=-\frac H6q_1+\frac{2v_2}3q_2-\frac{2H}3q_3,
\qquad F=-2q_2.
\]

逆包含もq1,q3がv2の倍であることから直接従い、**全 cokernel ideal は厳密に `(v2,F)`** となる。新しいGBは実行しない。

恒等式 `F=(v1+1)H−3` によりv1+1は商で単位である。`t=v1+1` とすれば

\[
v_1=t-1,\qquad v_0=-7t/4-1+3/(4t),\qquad v_2=0
\]

となる。従ってこの必要3jet parameter schemeはQQ上の滑らかなGmであり、被約・非空である。

## 4. 実際の有限 witness と CA20 の必要性

`v=(-2,0,0)` から低次18変数と高次18変数を独立に復元し、保存値と比較する。次にその具体的なB,C3と、第一 source `P0+aP1+a²P2` を**全く切り離した一つの有限5×5 pushforward**に代入する。a⁰,a¹,a²,a³ の全150出力が固定旗20空間に入ることを直接確認し、左逆で第二sourceの係数も復元して membership を検証する。

この witness は追加の CA binary 係数 index2をa²で−1だけ使う。同じ座標でその方向を欠くN19模型の不整合を、CA20へ転用することはできない。高次解の自由度によらないこの係数は

\[
(v_1+2)(3v_0+7v_1+5)/2
\]

と確認される。一方、全ての witness でこの係数が非零であるとは主張しない。

## 範囲と保存

受理するのは、第一sourceを既知の parabolic gaugeでfree10に取り、relative flag separationがa一次の ordinary chartでの、全二次解・全三次必要ideal・具体的な可解3jetである。可解な3jetから全次数 lift、代数的 pair branch、order4の結果、nonordinary pair の不存在、個別 conductor、B18原点の全 odd idealを推論しない。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/ambient_CA_origin_self_flag_3jet_review_20261103_v63.py \
  --output work/ca-self-flag3-independent-review
```

保存先は `reports/AMBIENT_CA_ORIGIN_SELF_FLAG_3JET_20261103_v63_REVIEW/report.json`。全 source / jet / witness の検証と入力hashを記録する。固定旗空間と parabolic gauge の長い既存証明は `AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW` を明示依存として保持する。
