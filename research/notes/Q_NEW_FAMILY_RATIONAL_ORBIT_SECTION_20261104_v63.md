# 五固定source座標上の有理Sp4 orbit代表

2026-11-04。u²=1600r, X=s/(48u), c=xi B u/40, W=(q+1)/3072と置く。
既存sectionのvは−32us/h、h=s²−25rq²となる。元の16次元表現の第一Cartan weightを
w1_iとし、Sp4のtorus g_u=diag(u,u^-1,1,1)を使うと

    beta_bar_i=u^(w1_i) beta_i

である。全16有理式をu²=1600rでexactにreduceすると、uの奇数部分は全て零となる。
従ってbeta_barは五つの有理座標(xi,r,s,q,B)上の**有理な16-vector**であり、十座標だけが
非零となる。有限torusの元なので全Sp4不変式を保つ。平方根を有理関数と仮定していない。

元のsourceと同じorbitを表すことからSp4-saturationは同じB_NEWに稠密である。
実target指数1の独立gateを使えば、これは一般orbitに対する有理代表の選択になる。
このpacket単独でgeneric stabilizerの自明性や十五次元hypersurface自体の有理性を
主張しない。全二十五次元Poisson componentの有理性とも区別する。

全係数は同名reportsのcertificate.jsonへ保存する。新しい点・高次数不変式・
source26の計算をせず、元の十Lie行列、全16座標と二次étale算術で恒等式を検証する。
source変数の名前sはここでは48uXを意味する。
