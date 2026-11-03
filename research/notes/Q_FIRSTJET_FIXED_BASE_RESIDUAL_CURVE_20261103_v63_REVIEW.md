# 固定 base の残余曲線と m 次数上界8：独立監査

2026-11-03。対象は `Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63`。
既に受理済みの二 adjoint 像・有限正規形・NEW25 の octic base 支配性を使う。
author code は import/再実行せず、新しい β 点、高次数点評価、GB は用いない。

## 1. 全係数と一点の微分

レビュー source は 4×4 行列式を24項の Leibniz 展開で計算し、固定 base の
全五 pencil 係数を照合する。p3 から a を決める分母は 2x(x²−y²)、p4=x²y²。
F4=20(x²+y²) と合わせ、x,y の選択は base の一般点上の高々8重の有限拡大。
その後の x,y,a,p1,p2 は五つの独立な base 座標になる。

p2 の式を u=gj、v=hi、m=u−v について解くと primary の二つの一次式になる。
g,h,u,v≠0 では同時符号の商は厳密に R=g/h、T=gh で、逆像は g²=RT の
二平方根である。レビューは二つの Q=RTL を g,h へ引き戻し、元の p1 と
determinant の式との差が全係数で零であることを確認する。平方根を有理関数と
扱っていない。

全斉次境界四辺と四隅を独立に照合する。二式は二重次数 (2,2)。Q2 は k1 に
依存せず、Q1 の k1 係数は −RT/2 である。

m の trace 式は graph によらない全対角 D chart（a,c,d,f を独立にする）で
確認する。同じ既存 QQ 点[0]だけで、4×4 Fraction 行列の forward dual 積に
より20偏微分を作り、既存20×16微分へ合成する。これは author の各項を
置換する微分とは別の実装である。六つの universal 不変量の minor は非零で
544/374 bit、全10軌道方向を殺す。NEW25 source へ合成した第六行は最初の
五行の線形結合であり、NEW25 の rank6 は主張しない。

## 2. 正確な次数評価

係数体を先に k1 を除いて取り、その後 k1 を超越元として加える。Q2 の各
既約因子は純超越拡大でも既約である。Laurent ring で Q1 と共通因子があれば、
Q1 の k1 係数である単元 −RT/2 も割ることになり矛盾する。従って generic
field 上で共通曲線成分はない。

xycfuv≠0 では Q1 の境界零点は四隅だけであり、そこでは Q2 の値
fu², −cf, −v², c がすべて単元。従って P1×P1 での共通点はすべて torus 上に
ある。proper な二 divisor の交叉長さは (2,2)·(2,2)=8。この数は総スキーム長さで
あり、八つの相異なる点や既約性を意味しない。

一つの実際の β 点で五 base 不変量と m の微分が独立だから、実際の既約像上で
その六つの有理関数は代数的独立である。従って実際の像の一般 base fiber 上で
m は非定数。有限 ordered-eigenvalue 拡大後の各選ばれた実成分でも非定数である。
c,f,u,v は係数体上の非零一次式なので、その成分上で恒等零とはならない。
u,v≠0 は g,h の非零も保証する。従って使用した Laurent chart はその一般成分に
交わり、m での generic fiber は上の長さ8の交叉の部分である。

よって、相対符号商と明示した有限 base 拡大の後、その選ばれた曲線の
m-line への一般次数は高々8である。全交点がその成分を占めること、正確な次数8、
全基礎体上の次数、NEW25 の base 上の multisection 次数は導かない。

## 3. 種数・spectral curve との区別

この証明は残余曲線の種数、有理性、楕円性を与えない。m 方向の compactification、
選択成分、正規化、分岐が必要である。一般の spectral plane quartic の genus や
その Prym を、そのままこの曲線に移すことはできない。特に通常の smooth
Sp4 spectral quartic の involution 商は genus1で Prym は次元2であり、今回の
一次元残余曲線をその Prym 自身と同一視する根拠にはならない。

この監査で受理するのは有限 recipe と次数上界だけである。NEW25 の未知の
追加 divisor、正確な multisection 次数、genus、Prym 同定は未解決のまま。

## 4. 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_firstjet_fixed_base_residual_curve_review_20261103_v63.py
```

測定時間・最大 RSS・全入力 hash・独立証明書は同名 REVIEW folder に保存する。
