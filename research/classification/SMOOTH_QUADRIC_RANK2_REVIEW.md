# 滑らかな二次曲面の階数二半単純cofactor：独立確認

2026-11-07。課題 `E_reduced_nonnormal`。
`SMOOTH_QUADRIC_CUBIC_SUPPORT.md` §6の二分岐を独立に確認した。

`code/classification_smooth_quadric_rank2_review.py` は、Gaussian rationalを二つのFractionで実装し、
QQ(i)上の疎多項式演算・消去を行う。producerとSymPyはimportしない。
`jobs/smooth_quadric_rank2_review.json` のgate後、約0.68秒でPASS。
全checker bytesと四入力reportのhashを
`reports/CLASSIFICATION_SMOOTH_QUADRIC_RANK2_REVIEW_20261107.json` に保持する。

q=x0²+x1²+x2²+x3x4、S=diag(0,0,0,1,−1)について、
全46不変係数に対する発散・q-normal条件を全200係数slotから組み立てた。
係数行列と拡大行列はともに階数35。
表示されたℓ三変数とtrace−1のA八変数は全条件を満たし、十一方向は独立なので全線形源である。

全strict Jacobiを再構成し、各係数がℓに一次、Aにaffine一次であることを確認した。
ℓ=0では五Pfaffian minor全てが恒等零。
SO3の非零非等方・等方の代表では、それぞれ全線形条件の階数7、拡大階数8となり解はない。
非等方代表ℓ=(a,0,0)のa≠0を1とできるのは、Jacobiがℓに厳密に一次だからである。
等方の非零ベクトルは複素SO3の単一軌道であり、零ベクトルは別に保持した。
従ってSだけをcofactorとするrank四支持は空である。

mixed型は零三次元blockのNを表示行列から独立に作り、
`Π_N=−5qN/6−E∧(Nx)/3` の発散・q-normal・N不変性を確認した。
全十一源へのN不変条件は係数・拡大とも階数8。
保存された代入を施した全三源(l1,m4,m7)は独立であり、全解を尽くす。
その全Jacobi係数のQQ(i)線形包は、厳密に

    3l1+5i,  l1(3m4+1)

の二式の線形包に一致する。
従って全Poisson源はl1=−5i/3、m4=−1/3、m7任意。
全m7についてJacobi恒等零と、実退化式

    D=(5/6)(x0+ix1)(x0²+x1²+x2²)q

およびその五つの偏微分とPfaffian minorの一致を確認した。
三次因子が必ず可約なので、滑らかなqと既約三次cの積という対象層から除外できる。

Poisson性から[S+N,Π]=0を得た後、線形表現のJordan分解により
[S,Π]=[N,Π]=0を課してよい点も整合する。
mixedの全一源自体は別因子型H+Q_rank3+Q_rank5であり、
このreviewではその成分帰属を主張しない。
