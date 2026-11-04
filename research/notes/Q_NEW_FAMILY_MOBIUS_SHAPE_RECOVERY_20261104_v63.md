# 二つのtarget比からBとqを有理的に復元する

2026-11-04。受理済みの短い七つの不変式写像から、**Bとqをtarget比と残るr,sで有理的に復元した**。既知τ-fixed source体を新しい五つの有理座標で表し直す結果であり、actual NEW target体との指数1や定義方程式は主張しない。

前packetの座標(ξ,r,s,q,B)と

\[
 h=s^2-25rq^2,\quad F=100r(r-1)q^2-4rs+s^2+4r,
 \quad F_4=(2/405)\xi^2F
\]

を使う。正規化targetを

\[
 Z_0=p_0/F_4^2,\qquad R=(p_3+4p_4)/p_0
\]

と取る。canonical恒等式からp₄/F₄²=1/1600−Z₀/4であり、p₃/F₄²=Z₀(R+1)−1/400。

圧縮された全p係数を直接組み合わせるとBが消え、

\[
 R=\frac{80qr^2-6r^2s+6r^2-11rs+10r+s}
 {8r(r-1)(5q-s+1)}.
\]

これはqについてMöbius変換で、その係数determinantは

\[
 -40r(r-1)^2J,\qquad J=(10r-1)s-10r,
\]

という非零多項式である。T,L,Nを

\[
 T=(r-1)R-2r,\quad L=8Rr(r-1)-6r^2-11r+1,
 \quad N=(s-1)L-(r-1)
\]

とすれば

\[
 \boxed{q=N/(40rT)},\qquad
 5q-s+1=(r-1)J/(8rT).
\]

さらに既存のBの復元式は

\[
 D=2359296000r^4(r-1)q^2(5q-s+1),\qquad
 \boxed{B=Z_0hF^2/D},
\]

で、q復元後はD=184320r(r−1)²N²J/T³となる。qとR、BとZ₀の**両方向の合成**が恒等となることをQQ有理式として検証した。従って表示したdense unit openで

\[
 \mathbb Q(\xi,r,s,q,B)=\mathbb Q(\xi,r,s,R,Z_0).
\]

sourceの未知shapeは(Z₀,R)を底としてr,sの二つに減る。

B⁰部分についてp₁,p₂,p₃,p₄,F₄²の定数係数線形関係を全係数で調べるとrank3で、kernel基底は(0,0,1/4,1,0)、(0,0,400,0,1)。p₁,p₂を含む新しい定数線形結合でBを消す方法はこの5列空間にはない。二次・有理な組合せを排除する結果ではない。

scopeは元sourceの非零分母に加えて、ξ,r,s,Z₀,T,r−1,N,J,h(q),F(q)を全て単元とする有理dense chart。これらの零集合への延長、写像全体のproperness、実際のtarget体の指数、残るr,sの有理復元は未証明である。

専用producerは`code/q_new_family_mobius_shape_recovery_20261104_v63.py`、保存は`reports/Q_NEW_FAMILY_MOBIUS_SHAPE_RECOVERY_20261104_v63/`。全小さい式、両方向の恒等式、定数部分rank、input hashを保持する。約0.34秒。新しい点、GB、resultant、source26展開は使わない。

先行の一般的な多変数rational cancellationは60秒で停止した。`code/q_new_family_two_shape_recovery_20261104_v63.py`と同名reportのexit142・frontierを保持する。そこで完了した復元式だけを本packetで独立に短く再現・保存し、失敗した約分部分を再試行していない。残る短い二shape式は別方式`Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63`で扱う。
