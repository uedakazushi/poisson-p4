# 下位の零一次部分層を除外する strict stabilizer 条件

2026-11-01。次の小さな障害は、前ノートの「射影 Poisson tensor の安定化群が
4次元以上」という必要条件を強化する。

T01 の一般点の四次元トーラスは、射影類だけでなく Poisson tensor 自体を保つ。
従って任意の T01 点 P において

    dim {A∈sl5 : L_A P=0 in H⁰(P4,Λ²T)} ≥4.

これは行列の kernel 次元の上半連続性から従う。P の非零スカラー倍を選んでも
kernel は変わらないので、射影化された Poisson 空間でも正しく定義される。

さらに、一般点の Cartan subalgebra を Gr(4,sl5) で閉包を取る。
Grassmannian の properness と、[a,a]=0、[a,P]=0 の閉条件により、
すべての T01 点は strict stabilizer 内に4次元の可換な Cartan-limit subalgebra を持つ。
これは安定化群の次元だけより強い必要条件である。

ただし P∈image(Λ²a→H⁰(Λ²T)) が極限でも成立するとはここでは主張しない。
この評価写像の階数は境界で下がり得る。射影化された評価写像には base locus があり、
像への所属が自動的に閉条件になるわけではない。
正しい閉じた incidence は可換性と [a,P]=0 である。

具体例として

    P=∂x∧∂y+∂z∧∂w+∂y∧(−xz∂z+xw∂w)

を取る。これは R0 の G=xzw に対応し、β=0、全134式を満たす。
退化因子は x0⁵ であり、完全分解条件にも反しない。

コードは全25個の線形座標変換を、150係数の Euler-wedge 商へ射影した。
strict 作用の階数は21であり、scalar identity を除いた stabilizer 次元は24−21=3。
strict stabilizer は二つの y-shear と z,w の逆向き再尺度化で張られる。
Poisson tensor 自身のスカラー方向も加えた行列の階数は同じ21なので、
射影 tensor の stabilizer 次元は25−21=4である。
従ってこの点は旧来の粗い条件を満たすが、strict 条件によって **T01 から除外される**。

証明書：`reports/TORIC_STRICT_STABILIZER_20261101.json`。
再現コード：`code/toric_strict_stabilizer_20261101.py`。
両階数の非零21次小行列と strict kernel の全行列を保存した。
残る9次元以下の候補層に strict kernel とその可換部分代数条件を課すことが次の課題である。
