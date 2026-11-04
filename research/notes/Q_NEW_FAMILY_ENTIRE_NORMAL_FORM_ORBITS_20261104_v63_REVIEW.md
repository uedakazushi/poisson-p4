# normal form全域のorbit同定：独立概念監査

2026-11-04。sourceとtargetがともにuv(32u+v)≠0にある場合に限る。
独立に確認したwhole-open旗とfull sigma証明を用い、generic field等式だけから
特殊parameterのorbitを推定しない。

Hは唯一の三重hyperplane、Q|Hのvertexはz、homogeneous modular Nの像は全u≠0で
span(ey,ez)。affine変換とtensor scalarに対しNtarget M=κ M Nsourceなのでこの像は保たれる。
従って一般affine旗は全域で正当。u=±8のcompatible空間の跳びはこの旗に影響しない。

η=κfに関する四つの低bracket係数式から(η²−1)a=(η²−1)b=0。
x/w blockの可逆性からa,bは同時に零にならずη²=1である。divergenceのy,h成分から
U=ηuとなり、全parameterでU=±uしかない。ここでaやcを割らず例外を落としていない。

同じ係数式からl±のpullbackがη=1では各lineの非零scalar倍、η=−1では交換となる。
三つのlineが現れるかどうか、logarithmic residueが有理数になるかどうかを仮定しない。
full projective同値はQも保ち、vertexからのconicのdual pairing比はこのline rescaleと交換に
不変である。det q2=−2048u(32u+v)、自己pairing32u/(32u+v)がdomain全体で単元なので
J=(v/u+16)²が全域で必要な不変量となる。

U=±uとJ一致が許す四つの候補は既受理tau/sigmaで全て実現できる。
そのmatrixの分母はsource/targetの同じdomain上で非零。従って必要十分性をQQの
rational normal forms間でも、Cを含む任意のchar0 field上でも受理する。
一般のcoarse k点がk-valued normal formを持つとはしない。

z=v/u+16を使うとparameter ringはk[u±1,z,(z²−256)^−1]で、二生成元はu,zの独立sign。
偶奇から固定ringはk[A±1,J,(J−256)^−1]、元ringは1,u,z,uzのfree rank4。
J=0では二つの幾何点と各二重のparameter fiber、それ以外では四点である。
parameter fiberの非被約性をPoisson componentの厚みと混同しない。

結論は表示されたnormal form族のPGL saturation内でのorbit同定と有限parameter商。
全C58の任意のsmooth-Q点の網羅や、全成分のgood PGL quotientの構築ではない。
