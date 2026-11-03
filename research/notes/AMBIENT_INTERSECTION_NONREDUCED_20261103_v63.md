# 全被約交叉 N の上にある実交叉の厚み：指定 Q での非被約性

2026年11月3日。このノートは二つの被約成分 CA,CB の **実際の scheme-theoretic intersection**
についての主張である。正規化同士の tensor product を target の交叉と同一視しない。
指定 Q は継続している v24 a=1 の点であり、別の共鳴点a=2/9を対象にしていない。

## 定理

実交叉 D=CA∩CB は指定 Q で非被約である。
同じ結論は、7次元臨界 rank2 一次ジェット軌道の上の、共通の滑らかな10次元非正規層の全点で成り立つ。

具体的に S61 で

    dim T_Q(CA)=dim T_Q(CB)=28,
    dim T_Q(N)=26,
    dim T_Q(D)≥27

である。既に全被約交叉が N であることを証明しているので、最後の厳密な接空間の差が
交叉の非被約な厚みを検出する。
nilradical の全生成元、冪零次数、全 primary decomposition はまだ決めていない。

## 1. 正規化の微分だけを使ってはいけない

受理対象の全完成局所環表示は、Nが A7、CA,CB が A8、ただし

    Ac=C[[u10,s1,…,sc,p²,p³,s1p,…,scp]]

というものである。
その正規化の通常の一次微分は u10,s1,…,sc の d−1方向しか見ない。
target の真の cotangent は、これらに加え p²,p³,sip を独立な極大 ideal 生成元として持つ。
従ってCA/CBの埋め込まれた28次元 tangentは、正規化からの17/18次元微分像ではない。

## 2. 共通の七次元群因子を外した全 cotangent

群作用で b0=(τ,A,B,C,D,F)=(0,0,1,1,0,0) へ移す。
共通の七次元 smooth orbit factor を外すと、全非正規点は

    p=s=r=k1=k3=k4=k6=k7=k8=0,
    (k0,k2,k5)=(c0,c2,c5)

で表せる。Nではもとよりr=0である。
実際の全61出力は保存済みで、c0,c2,c5 はこの計算では独立な記号として扱う。

target から選んだ regular 座標 K と s を固定し、正規化の source fiber を

    w(p)=w0+w1p+w2p²+w3p³ mod p⁴

と展開する。選択出力の係数行列を L(p)、定数項を b(p) とすると、
L0 は可逆なので

    wn=−L0^−1*(bn+Σ_(i=1)^n Li*w_(n−i))

で全補正係数が正確に決まる。
これを全61出力へ代入したp²,p³係数が、対応する二つの cotangent 生成元に対する埋込み tangent vector を与える。
さらに K,s の通常微分と、K*p,s*p の mixed coefficients を加える。
これらは上の Ac の極大 ideal の全生成元を尽くすので、全 target tangent space を張る。
正規化上の適当な曲線の数本から下界だけを取ったのではない。

## 3. QQ(c0,c2,c5) 上の正確な行列

この方法で61行の行列を三つ作った。七次元 orbit factor を除いた階数は

    rank TN=19,
    rank TA=rank TB=21,
    rank [TA TB]=22

となる。TN はTA,TB双方に含まれる。
全行列の係数は c0,c2,c5 の次数1以下の多項式である。
QQ(c0,c2,c5) 上で union rank22 であることは、全23次小行列が多項式として恒等的に零であることを意味する。
従って全ての特殊化で union rank≤22となる。

各成分の階数21とNの階数19は、その全局所環Acによって全ての非正規点で一定である。
よって各特殊化で

    dim(TA∩TB)≥21+21−22=20>19=dim TN.

共通の七つの smooth orbit directions を戻すと、定理の28,28,26,≥27が得られる。
generic な三係数の標本だけから指定 Q に結論を移していない。
上界が全パラメータで成立する多項式恒等式であることを用いている。

## 4. 実 target 交叉の非被約性

閉部分スキームの scheme-theoretic intersection の tangent space は、
同じ ambient tangent space 内の二つの tangent spaces の交叉に厳密に等しい。
ここでは同じS61の全61座標内でTA,TBを作っているので、上の交叉は実際のDの接空間である。

全被約交叉の定理により Dred=Nである。
もし D が対象の点の近傍で被約ならば、その局所環はNの局所環に一致し、
接空間も同じ26次元になるはずである。
しかし実際は少なくとも27次元なので矛盾する。従って D は非被約である。

この推論はB_A tensor B_Bの余分な冪零元をtargetへ降ろしたものではない。
個々の被約target環の **全cotangent** をそれぞれ実ambientに埋め込み、その和idealの一次項を直接比較している。
全intersection idealの消去計算は今後も必要である。

## 証明書

* `code/ambient_intersection_tangent_excess_20261103_v63.py`
* `reports/AMBIENT_INTERSECTION_TANGENT_EXCESS_20261103_v63/report.json`

全記号行列、三つのsource点、rank21/21/19、union rank22、包含は QQ(c0,c2,c5) 上で検査し、
三つのQQ特殊化も照合した。測定1.67秒。
全被約交叉・全成分局所環の二つの入力定理と共に監査して使う。
