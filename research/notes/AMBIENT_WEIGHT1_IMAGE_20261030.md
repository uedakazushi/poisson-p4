# 二つの既知28次元成分の共通一次jet像

2026年10月30日。対象は `legacy/v06/ambient_families.tex` の既知成分
\(\mathcal C_A,\mathcal C_B\)。新しい成分の主張ではない。

## 結論

両成分に対応する S61 成分の、原始 gauge における一次jet射影の像は、
同じ9次元の閉既約部分集合 \(Z\subset V_{16}\) である。特に両者の像は係数写像
\(K_\beta:U\to\Lambda^2_0U^*\) の階数3以下の部分に入る。
ここで像は複素点の集合とその被約閉部分集合を意味し、スキーム論的な像の構造は主張しない。
一般の一次jetに対する [T01 の優越性](TORIC_WEIGHT1_DOMINANCE_20261030.md) と対照して、
CA/CB は一般の一次jet上には現れない。

既存の28次元という成分次元と S61 の次元差9を用いると、両者の19次元の S61 成分の
射影 \(\mathcal C_{A,S61},\mathcal C_{B,S61}\to Z\) の一般ファイバー次元は10となる。
これはファイバーの既約性・同型類や、二成分の交叉次元を決める主張ではない。

## 共通の六次元線形族

座標を \((x,y,z,w)\)、定数 symplectic bivector を
\(J=\partial_x\wedge\partial_y+\partial_z\wedge\partial_w\) とする。
\(h_2=A xz+B xw+Cz^2+Dzw+Fw^2\)、
\(X_{h_2}=-(h_2)_w\partial_z+(h_2)_z\partial_w\) と置く。

両族の原始一次項は、パラメータを取り替えると共通の形

\[
L_0=t x(\partial_x\wedge\partial_y-\partial_z\wedge\partial_w)
 +\partial_y\wedge\left(X_{h_2}+\frac t2(z\partial_z+w\partial_w)\right)
\]

になる。CA の \(f(1,x)=1+px+O(x^2)\) では \(p=5t/2\)、CB では \(p=-5t\)。
したがって、二形式 \(\beta=-JL_0J\) の六次元線形空間 \(S\) は完全に同じである。
具体的には

\[
\beta_{12}=tx,\quad\beta_{34}=-tx,\quad\beta_{23}=\beta_{24}=0,
\]

\[
\beta_{13}=(h_2)_z+\frac t2w,\qquad
\beta_{14}=(h_2)_w-\frac t2z.
\]

これは各パラメータについて厳密な等式であり、標本の一致ではない。
任意の \(t,h_2\) は、原族で \(f=x_0^2+p x_0x\)、\(h=x_0h_2\) と選べば実現する。
評価点の値は常に J なので、rank4 の開条件を失わない。

## 任意のrank4点と定数項の正規化

両族の実際の退化五次式はそれぞれ \(x_0^3f\)、\(x_0f^2\)。したがって
rank4 点は必ず \(x_0\ne0\) にある。affine translation により任意のrank4点を原点へ移せる。
CA では族の形がそのまま保たれる。CB の Euler 場の平行移動による定数ベクトル項は
Hamiltonian 項へ吸収でき、三次以下の \(h\) の許容範囲を保つ。

原点での値を

\[
\Lambda=\begin{pmatrix}
0&a&0&0\\-a&0&d&e\\0&-d&0&g\\0&-e&-g&0
\end{pmatrix},\qquad ag\ne0
\]

とする。CA では \(g=1\)、CB では \(g=a\)。

\[
R=\begin{pmatrix}
1&0&0&0\\0&1/a&0&0\\d/a&0&1&0\\e/(ag)&0&0&1/g
\end{pmatrix}
\]

なら \(R\Lambda R^t=J\) である。この変換は x 座標と y への非依存性を保つ。
正規化後の一次項は

\[
L_{12}=sx,\quad L_{34}=rx,\quad L_{13}=L_{14}=0,
\quad L_{23}=Y^z(x,z,w),\quad L_{24}=Y^w(x,z,w),
\]

ただし \(Y\) は線形で \(\operatorname{div}_{zw}Y=-r\) を満たす。
ここで CA は \(r=0\)、CB は \(r=s\)。これは Jacobi の定数係数から従う。
\(\tau=(s+r)x\) による原始正規化後には

\[
t=(2s-3r)/5,\qquad (L_0)_{12}=tx,\quad (L_0)_{34}=-tx,
\quad \operatorname{div}_{zw}Y_0=t
\]

となる。従って \(Y_0-\frac t2 E_{zw}\) は divergence-free で、上記の \(X_{h_2}\) と書ける。
定数交差項 \(d,e\) を零と仮定する必要はない。この一般定数正規化もコードで記号的に検証した。

異なる射影座標系で同じrank4点を正規化した差は、その点の放物部分群から生じる。
perspective は原始 gauge で除かれ、残りは CSp4 の作用である。六次元空間 S は線形で、
係数のスカラー倍を許すので、その Sp4 飽和は CSp4 飽和とも一致する。

## 座標に依存しない閉部分集合

定数 symplectic 二形式を J の行列で表す。各非零ベクトル \(v\in U\) に対し
\(\ell_v=\iota_vJ\in U^*\) とし、

\[
S_{[v]}=\{\beta\in V_{16}: K_\beta(v)=0,\quad
\iota_v\beta(y)\in\mathbb C\,\ell_v(y)\ell_v\}
\]

と定める。第二条件は、y について線形な一形式全体の空間で、表示した一次元部分空間へ
入るという条件である。定数倍で v を替えても条件は変わらない。

\(v=\partial_y\) の場合、これらの条件は
\(\beta_{12}=tx\)、\(\beta_{34}=-tx\)、\(\beta_{23}=\beta_{24}=0\)、
\(\beta_{13},\beta_{14}\) が x,z,w の線形形式であること、および閉性の一関係を与える。
したがって \(S_{[\partial_y]}=S\)、次元は6。
コードでは24係数上の制約行列の階数18と、S の六つの基底によるその核の生成を検証した。

Sp4 は \(\mathbb P(U)=\mathbb P^3\) に推移的に作用する。よって

\[
\mathcal I=\{([v],\beta):\beta\in S_{[v]}\}\subset\mathbb P^3\times V_{16}
\]

は階数6のベクトル束の全空間で、既約・9次元である。条件は各固定 v に対して線形で、
\(\ell_v\ne0\) なので相対的な階数も一定。第二条件は二テンソルの比例性という閉条件であり、
\(\mathcal I\) は表示した積の中で閉である。

射影 \(\mathbb P^3\times V_{16}\to V_{16}\) は proper なので
\(Z=\operatorname{pr}_2(\mathcal I)\) は閉・既約で、次元は9以下。
\(K_\beta(v)=0\) より常に \(\operatorname{rank}K_\beta\le3\)。

有理点 \((t,A,B,C,D,F)=(2,3,5,7,11,13)\) では係数写像の階数は3で、核は
\(\mathbb C\partial_y\)。このため一般の像点の v は一意であり、\(\dim Z=9\)。
さらに簡単な証人は
\(\beta_*=x\,dx\wedge dy+(z+w)dx\wedge dz+w\,dx\wedge dw-x\,dz\wedge dw\)。
\((t,A,B,C,D,F)=(1,0,0,1/2,1/2,1/2)\) に対応し、係数写像の
行 \((12,13,14)\)、列 \((x,z,w)\) の小行列式が1である。
別の証明書として、六パラメータ微分と十個の Sp4 作用を合わせた行列にも
非零の9×9小行列式 \(-2016\) がある。
上界はこの有限個の階数標本から推測したものではなく、上の incidence から得ている。

Sp4 飽和 \(\mathrm{Sp}_4\cdot S\) は Z と等しく、S の全点が両原族で実現するので、
Z の全点は両成分から実現する。逆に原族の全rank4点とその射影座標変換から得る
原始jetは Z に入る。Z は閉なので、対応する S61 成分の閉包を取ってもこの包含は保たれる。
従って両成分の射影像は同じ Z である。

## 範囲・再現

この結論は、CA と CB が同じ成分であることや、同じ一次jetの二つの構造が等しいことを
意味しない。成分の相異性は従来の退化因子の重複度による。一次jetだけでは両者を分離できず、
その区別は高次係数のファイバーに残る。

`code/ambient_weight1_image.py`、`reports/AMBIENT_WEIGHT1_IMAGE_20261030.json` に、
六次元族の一致、一般定数正規化、intrinsic fiber の階数、放物部分群の作用、
9次元像・3階数係数写像の厳密小行列式を保存した。

再現：`.venv/bin/python code/ambient_weight1_image.py --output reports/<new-name>.json`。
既知の成分数の下界57は変更しない。
