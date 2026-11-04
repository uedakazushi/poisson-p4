# 追跡因子による二shape有理式：直接implicitizationの有限な入力

2026-11-04。新しい次数の不変式・点を作らず、受理済みの七つの短い写像へqとBの復元式を代入した。一般的な多変数gcd/cancelを使う最初の方法は60秒で停止したが、**qの多項式斉次化と既知因子だけの完全除算**により、残る二shapeの三つの局所化方程式を約0.81秒で構成できた。implicitization自体は未実行である。

## 座標と局所化

`Q_NEW_FAMILY_MOBIUS_SHAPE_RECOVERY`の記号を使う。

\[
 T=(r-1)R-2r,\quad N=(s-1)(8Rr(r-1)-6r^2-11r+1)-(r-1),
 \quad J=(10r-1)s-10r.
\]

q=N/(40rT)、B=Z₀hF²/D。さらに

\[
 H_n=64rT^2s^2-N^2,\quad
 F_n=(r-1)N^2+16rT^2(s^2-4rs+4r),\quad\alpha=2/405
\]

とすれば、直接の恒等式として

\[
 h(q)=H_n/(64rT^2),\quad
 \alpha F(q)=\alpha F_n/(16rT^2),
\]

\[
 L_0(q)=\alpha^2\,184320r(r-1)^2N^2J/T^3
\]

を得る。L₀はp₀のξ⁴B/h係数である。

## 正確な代入方法と短い係数関数

既存のP(r,s,q)に対しm=deg_qPを取り、

\[
 \widehat P=(40rT)^m P(r,s,N/(40rT))
\]

を項ごとの多項式乗算で構成した。r,T,r−1,N,J,Hₙ,Fₙという**名前を固定した因子だけ**で完全除算し、除算ごとに商×因子が元へ戻ることを確認した。多変数gcd、一般的なcancel、factorization、GB、resultantは使わない。係数のscalar contentだけは正規化する。

正規化targetをZ₁=p₁/F₄²、Z₂=p₂/F₄²、Z_N=N₁₂/F₄³とすると、結果は

\[
 Z_1=A_1+Z_0C_1,\qquad Z_2=A_2+Z_0C_2,\qquad
 Z_N=Z_0E_1+Z_0^2E_2
\]

である。各分子は(r,s,R)のprimitive整数多項式で、定数倍を含めてcertificateに全係数を保存した。

|係数関数|分子項数|分子の(r,s,R)多重次数|定数倍|分母|
|---|---:|---|---:|---|
|A₁|33|(5,2,2)|1/50|Fₙ|
|C₁|15|(4,1,1)|1/8|r(r−1)²J|
|A₂|212|(10,4,4)|1/400|Fₙ²|
|C₂|17|(4,1,1)|1/16|r(r−1)²J|
|E₁|79|(6,3,3)|1/20|FₙJ|
|E₂|41|(5,2,2)|5/4|r(r−1)²J²|

たとえばA₁は33項の保存分子を50Fₙで割ったもの。他も同じ規約である。元のq由来のN²、Hₙ因子は全て明示的な完全除算で消えた。

## 等価な二shape方程式と次の計算費用

名前付き分母の最小公倍積を掛け、Z₁,Z₂,Z_Nを独立なtarget変数として得た方程式は次の大きさとなった。

|方程式|shape次数(deg_r,deg_s)|R次数|target・Z₀を別変数と数えた項数|
|---|---|---:|---:|
|Z₁−A₁−Z₀C₁=0|(9,3)|3|284|
|Z₂−A₂−Z₀C₂=0|(14,5)|5|1076|
|Z_N−Z₀E₁−Z₀²E₂=0|(10,4)|4|471|

対応する分母はそれぞれr(r−1)²JFₙ、r(r−1)²JFₙ²、r(r−1)²J²Fₙ。
これらは**追跡済み因子を約分した局所化等価式**であり、未調査の因子まで含むglobal primitive/gcd最小形とは主張しない。

対象openは

\[
 \xi r s Z_0(r-1)T N J H_nF_n\ne0.
\]

全ての除外因子を保持し、零集合での延長や全chartのpropernessを仮定しない。既存の特性基底への優越と、Z₃=p₃/F₄²=Z₀(R+1)−1/400から、Z₀,R,Z₁,Z₂は代数的に独立。従ってQ(Z₀,R)上の二変数(r,s)から(Z₁,Z₂)への**一般的な有理写像は有限な函数体拡大**を与える。このことだけでは次数・ramification・r,sの有理回収は決まらない。

最初の二式はsについて3次と5次。将来の一つのSylvester resultantは8×8で、r次数上界は9·5+14·3=87、R次数の素朴な上界は30である。ただし境界因子・base locusと局所化を保持した評価が必要。resultant、pseudo-remainder、quotient algebra、第三式による消去は本packetでは実行していない。先にA₂とA₁²のFₙ上のpoleを比較して第二式を小さくできる可能性があり、これも別の小さい検査とする。

将来得られる低不変式間のrelationがNEW25で消えるとしても、普遍的に独立な六つの低生成元からのpullbackがNEWの既約定義式そのものであるかは、次数一致・因数分解・対象成分の認定が別途必要。actual target体の指数1を仮定しない。

## 保存・失敗履歴・再現

成功した方法は`code/q_new_family_two_shape_homogenization_20261104_v63.py`と
`reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/`。全homogenized係数、完全除算回数、係数関数、三つの局所化式、frontier、入力hashを約0.51MBのcertificateに保存した。60秒/2GB/32MiB cap内、一回約0.806秒。

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_two_shape_homogenization_20261104_v63.py \
  --output work/two-shape-homogenization-replay
```

別の新しいoutputを指定する。先行するgeneric rational cancellationの60秒exit142は、`Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63/report.json`とそのcode/frontierに失敗として残した。成功済みq/B回収は専用Mobius packetで短く独立再現した。失敗した方法を同じ条件で再試行せず、source26・新点・高次数の展開も行っていない。
