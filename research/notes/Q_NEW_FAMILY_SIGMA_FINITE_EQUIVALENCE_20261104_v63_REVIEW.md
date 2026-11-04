# 実際のsigma有限変換とfull PGL moduli体の独立確定

2026-11-04。対象は凍結 `Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63`。
producerのaffine SymPy展開・cancelとは別に、QQ(u,v)を係数体とする
五変数(h,x,y,z,w)の疎多項式環で、元の全affine Pを直接homogenizeする。
保存された二つのA,tからMを組み、各epsilonについて全150 homogeneous係数を照合する。
producerのimport・再実行、同値solver、点やprimeの評価は行わない。

p=32u+vとすると、双方について独立な置換展開により

    kappa M rawP(u,1,v)(X) M^t = rawP(u,1,−p)(M X),
    det(M)=det(A)=−v^7/p^7,  kappa=−p²/v²

がQQ(u,v)上で成り立つ。行列のdeterminantも順列和で計算した。
分母はu,v,pと非零有理数だけであり、表示された `u*v*p!=0` 全体で可逆である。
従って単なる一般点のsampleではなく、このopen全体のactual full projective Poisson
同値である。raw equalityはEuler quotientへ下降し、kappaは射影化で消える。
sigma後にはpが−vとなるので、このopenはsigmaでも保たれる。

以前の独立下側gateは、評価点とscalarにも依存しないfull PGL不変量
`A=u²,J=(v/u+16)²` を構成している。その凍結reviewを引用し、こちらでは
変更しない。既知tauと今回のactual sigmaを合わせると

    k(A,J) ⊂ E_full ⊂ k(A,b),  b=uv,
    sigma(b)=−32A−b.

右端におけるsigma固定体はちょうどk(A,J)。従って

    k(X_NEW25)^PGL5 = k(u²,(v/u+16)²),  k=QQ,C.

これはfull Poisson成分のmoduli体の**明示的な純有理二変数表示**である。
firstjet像B15の五次元Sp4商とは異なる。

幾何的generic parameterでA,Jの同じ値を持つ組は
`(u,v),(-u,-v),(u,−32u−v),(-u,32u+v)` の四つだけであり、全てactual変換で
実現される。tauとsigmaはparameter上で可換なinvolutionなので、一般parameter間の
同定は正確にKlein fourである。この結論は一つのtensorの有限isotropy群の計算、
特殊parameter・境界の網羅、全不変式環の生成、定義方程式の構成を意味しない。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_sigma_finite_equivalence_review_20261104_v63.py`
（一回60秒/2GB上限）。
