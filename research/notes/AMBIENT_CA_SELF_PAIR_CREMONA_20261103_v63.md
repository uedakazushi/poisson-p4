# CA の実 self-pair 族：一様な Poisson Cremona・有理葉・既知成分への所属

2026-11-03。`AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63` の実際の多項式族 P(a)を対象とする。a≠0の全点を、未知のparameter因子を反転せずに扱う。新しい成分を主張するものではない。

## 1. x=1 の一様な多項式 Poisson 座標

全150係数で保存された canonical tensorは

\[
P_a=\partial_y\wedge V_a,
\quad V_a=(x^2+12ahx/5)\partial_h-3ax^2\partial_x/5
 +(hx+2axz/5)\partial_z+(h^2+z^2-8awx/5)\partial_w.
\]

x=1の射影chartで yもy/xの意味にすると、実際のPoisson tensorは

\[
\partial_y\wedge D_a,\qquad
D_a=(1+3ah)\partial_h+(h+az)\partial_z+(h^2+z^2-aw)\partial_w
\]

である。Vのx成分によるprojective補正を引いている。

a≠0で

\[
H=h+\frac1{3a},\qquad Z=z-\frac h{2a}-\frac1{2a^2}
\]

と置くと、`D(H)=3aH`、`D(Z)=aZ` となる。次の**二次多項式**を取る。

\[
\begin{aligned}
Q_a(H,Z)={}&\frac{4a^2+1}{28a^3}H^2+\frac1{5a^2}HZ+\frac1{3a}Z^2\\
 &+\frac{1-2a^2}{12a^4}H+\frac1{3a^3}Z+\frac{a^2+1}{9a^5}.
\end{aligned}
\]

`h=H−1/(3a)`、`z=Z+H/(2a)+1/(3a²)` を代入すれば厳密に

\[
\bigl(a(3H\partial_H+Z\partial_Z)+a\bigr)Q_a=h^2+z^2.
\]

二次以下の各単項式の固有値は、aのそれぞれ7,5,3,4,2,1倍であり、a≠0では共鳴しない。従って

\[
W=w-Q_a(H,Z)
\]

とすれば `D(W)=−aW` となる。これはJacobian determinant1を持つ二次の多項式自己同型で、逆は上のh,zと `w=W+Q_a(H,Z)`。全Poisson tensorが正確に

\[
\partial_y\wedge(3aH\partial_H+aZ\partial_Z-aW\partial_W)
\]

へ移る。さらにY=y/aとすればaも消せるが、これは有理Poisson同値であってPGL同値の主張ではない。

## 2. P4 上の実際の二次 Cremona

\[
A=h+\frac{x}{3a},\quad
B=z-\frac h{2a}-\frac{x}{2a^2},\quad
C=xw-x^2Q_a(A/x,B/x)
\]

とする。Cは同次二次式である。上の自己同型は射影的に

\[
[h:x:y:z:w]\dashrightarrow[Ax:x^2:xy:Bx:C]
\]

と書け、逆写像は

\[
[H:X:Y:Z:W]\dashrightarrow
[X(H-X/(3a)):X^2:XY:
 X(Z+H/(2a)+X/(3a^2)):
 XW+X^2Q_a(H/X,Z/X)].
\]

両合成はそれぞれx³とX³倍の恒等写像となる。Cのz²係数は−1/(3a)なので、全座標に共通のx因子はなく、全a≠0で次数2である。

同次次数1の有理表示 `(A,x,y,B,C/x)` はEuler場を保ち、全5×5 Jacobianによるactual tensor pushforwardも確認した。その像は、上のaffine diagonal tensorを斉次化してcanonical化した二次tensorに等しい。葉の分布だけの同値や、未知の函数倍への置換ではない。

## 3. 完全な有理 Casimir 体と一般葉の閉包次数

二つの有理Casimirは

\[
f=H/Z^3=Ax^2/B^3,\qquad g=ZW=BC/x^3.
\]

`c1=f,c2=g,t=Z,Y=y` とすると逆は

\[
\begin{aligned}
h&=c_1t^3-1/(3a),\\
z&=t+c_1t^3/(2a)+1/(3a^2),\\
w&=c_2/t+Q_a(c_1t^3,t),\qquad y=Y.
\end{aligned}
\]

全函数体は `C(c1,c2,t,Y)` となり、Poisson tensorは `at∂Y∧∂t` である。従って有理Poisson中心は**正確に** `C(f,g)` である。c1c2≠0の一般葉は、このchartで `C×Gm` となり、そのZariski閉包は有理曲面である。

y方向を射影する商P3の曲線を、上の四座標 `(h,1,z,w)` にtを掛けて表示すると

\[
[c_1t^4-t/(3a):t:
t^2+c_1t^4/(2a)+t/(3a^2):
c_2+tQ_a(c_1t^3,t)].
\]

最後の座標のt⁷係数は `c1²(4a²+1)/(28a³)`、t⁶係数は零、t⁵係数は `c1/(5a²)` である。c1c2≠0では有限t≠0で第2座標が非零、t=0で最後の座標c2が非零、無限遠でも最高次係数が非零である。従って同次化にはbasepointがなく、t=B/xを戻せるので像にbirationalである。

よって一般葉のP4閉包はvertex `[e_y]` を持つ有理曲面coneで、次数は

* a≠0、4a²+1≠0では **7**。
* 4a²+1=0の二つの非零複素parameterでは **5**。

特殊なCasimir値c1c2=0の葉をこの一般次数で記述していない。

## 4. LPB39 と Log(1,1,1,2)28 への actual 所属

**LPB:** 定vector factor∂yと、yに依存しない二次場Vがあるので、旧 `legacy/v14/code/rank2_lower_bounds/lpb_family.py` の定義に直接入る。旧中心e4へyを送る置換 `(0,1,4,2,3)` を用い、raw及びcanonical全150係数で `e4∧V` と一致させた。従って全a、a=0も含め、既存LPB39に入る。

**Log1112:** a≠0では、divisorsを `(A,B,x,C)`、degreeを `(1,1,1,2)` とする。`dlog f,dlog g` のresidueはそれぞれ

```
(1,-3,2,0), (0,1,-3,1).
```

両者はEuler-horizontalである。これらのwedge係数は、旧 `verify_rank2_components.py` のlog builderでparameters `(1/4,-3/4,7/4)` を取ったものに正確に一致する。積を払った2formをω、`Euler∧P_a` のvolume-dualをΩPとすると、全10係数で

\[
\omega=-\Omega_P/a
\]

となる。旧builderの `dω` のdual tensorも全150係数で

\[
P_{\log}=-5P_a/a
\]

と一致する。したがって単に同じfoliationを持つという判断ではなく、actual tensorの非零定数倍として旧Log1112族に所属する。divisorsが特殊な場合も、旧sourceの多項式表示の特殊化はその閉包に入る。

既存Log1112が全Poisson空間の実28次元成分、LPBが実39次元成分であることは、旧 `legacy/v14/data/six_rank2_components.json` の受理済みtangent/family certificateを再利用する。その大計算は再実行していない。

P(a)はa=0でも正則・非零の多項式族であり、Log1112は閉だから、**a=0を含む全多項式直線もLog1112に入る**。元packetによるCA28所属と合わせ、少なくともCA28・LPB39・Log1112_28という既知三成分への包含がある。CBや他の既知成分への所属・非所属、全conductor、全self-pair locusをここから推論しない。

## 再現と保存

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/ambient_CA_self_pair_cremona_20261103_v63.py \
  --output work/ca-self-pair-cremona-replay
```

`reports/AMBIENT_CA_SELF_PAIR_CREMONA_20261103_v63/` に、Q・双方向Cremona・全Casimir逆式・全log係数・旧sourceの対応とhashを保存する。GBや数値補間は使わない。
