# 正 weight6 軌道に沿う実八変数形式模型と二次 obstruction

2026-11-03。正軸 eta1 の物理的 canonical tensor を P とする。
全126 canonical 方向と全350 strict Jacobi 係数を使い、23変数の滑らかな orbit 因子を除いた
**埋め込み次元8の実形式模型**を得た。選択95式の形式逆写像と、残る255式の残余形式級数による
全次数の定義を持つ。有限の八変数多項式 ideal と同定した意味ではない。

その二次 obstruction 空間の次元は正確に5である。二次式だけの ideal の radical は二つの5平面の和集合、
次元5・次数3である。従って実形式横断環の次元は5以下、全 projective Poisson 局所次元は27以下となる。
実次元が5であることや、この二平面が実枝へ持ち上がることは主張しない。
全 eta≠0 点は一つの既証 projective 軌道なので、同じ局所結論が軸の全非零点に成り立つ。

## 全方程式からの形式積

canonical 線形空間の基底 C は凍結済み全126 Smith 証明書を用い、C の独立 raw 座標行は恒等行列となる。
eta1 の P と、GL25作用および scalar の全150方向をその基底で表した。
実 orbit 行列の rank23から23本と23 raw 座標を選び、その minor

    −332306998946228968225951765070086144/14348907

が非零である。選択座標を中心値に固定した103次元 affine 線形切断は実 orbit に横断する。
選択 Lie 方向の局所作用と scalar による写像の微分は同型なので、形式逆写像定理により
全 canonical Poisson 完備局所環は切断環に23形式変数を加えた環になる。
射影化で非零 scalar 一本を除くと22個の滑らかな因子となる。
群座標の大域一意性や別の点で同じ chart を使えるとの主張は必要ない。

この103方向において、全350 Jacobi 線形化の rank は95である。
95列・95行の厳密 QQ 逆行列 B^(-1) を保存し、B*B^(-1)=I95を検証した。
残る八座標を z0,...,z7 とする。線形 kernel T とその物理的150×8像を保存し、
全350微分が消え、23 orbit cut も消えることを確認した。
したがって canonical affine 接空間次元31、projective 接空間次元30、形式横断の埋め込み次元8である。

95選択式は dependent95変数を八自由座標の形式級数として一意に解く。
その形式解を**残るすべての255元式**へ代入することで、実八変数形式 ideal が定まる。
収束・形式逆写像の一般定理による全次数の presentation であり、低次 jet を全次数模型とみなしていない。
canonical divergence 零の模型では strict Jacobi350と射影 Poisson 条件は既存 BV 正規化により同値である。

## 厳密な二次補正

物理的八方向 Tz を元の二次 Jacobi tensor に直接入れ、全350×36二次係数 Q を構成した。
95選択式の二次解は −B^(-1)*Q_selected であるため、残余二次係数は

    O=Q−J_dependent*B^(-1)*Q_selected

である。全選択95行の O が零になることも確認した。
厳密 QQ row reduction による五式は

    z1²+120z5z7−14z7²,
    (15z1z3−5z1z6+12z3z7+8z4z5−2z4z7+16z5z6−8z6z7)/15,
    (3z1z4+2z1z6+36z3z7−16z4z5+12z4z7−32z5z6+8z6z7)/3,
    (4z1z5+36z5z7−z7²)/4,
    z7(z1−20z5+5z7).

この二次 ideal I2 の小さな QQ Gröbner 証明書は次元5、次数3を与える。
K=(z1,z7,z5(z4+2z6)) とすると、全五式は K に入り、
z1³、z7³、(z5(z4+2z6))² が I2 に入ることを厳密に検査した。
Kは二つの相異なる線形素 ideal

    (z1,z7,z5) ∩ (z1,z7,z4+2z6)

の共通部分である。従って sqrt(I2)=K であり、その二つの被約5平面を得る。
実 ordinary initial ideal は I2 を含むので、実形式環の次元上界5は従うが、
それ以上の低次・高次 cancellation が新しい関係を生む可能性は残る。

## 保存と次の課題

`code/q_positive_weight6_quadratic_model_20261103_v63.py` は元の物理的 tensor の Jacobi を直接微分し、
全350式を保持する。`reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63/` に
全95逆行列、全350×103線形化、全150×8kernel、全350×36の未補正・補正二次係数、
および小さな `quadrics.m2`・実行ログを保存する。

最初の正確な QQ 計算は4.5秒だった。大規模 GB や blind な高次 jet を繰り返したものではない。
次の課題はこの八変数形式 presentation の真の三次・四次 initial を、低次 cancellation を含めて求めることである。
