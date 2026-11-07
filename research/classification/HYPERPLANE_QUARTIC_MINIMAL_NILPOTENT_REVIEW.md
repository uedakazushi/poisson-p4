# 最小冪零・translation 零の四次因子：全支持被覆の独立確認

2026-11-07。課題 `E_reduced_nonnormal`、純冪零 cofactor 表の型8。
`HYPERPLANE_QUARTIC_MINIMAL_NILPOTENT_SUPPORT.md` の全支持証明を独立に精査した。
結論は、単純 H と既約四次 F を持つ `Z=y∂x` 型の全源が既存 MINIMAL_M 成分へ入ることである。
他の Jordan 型、translation 非零型、全成分の尽くしはここから従わない。

## 完全性の数学的確認

`[Z,P]=0` と全90係数の源から、全括弧を

    Pxy=a, Pxz=xl+b, Pxw=xm+d,
    Pyz=yl, Pyw=ym, Pzw=f

と置ける。`a,b,d,f` は x に依存しない二次以下、`l,m` は一次以下である。
`F=Pf(P)=af+y(dl−bm)` は x に依存しない。
既約四次なら `F|y=0≠0`。この制限が零なら、四次斉次化も y で割れてしまう。
modular 恒等式 `P grad F=−FZ` の z,w 成分を y=0 に制限すると、
非零多項式 `f0` と `∂zF0,∂wF0` の積が零である。
従って F0 は非零定数、積 `a0f0=F0` の両因子も多項式環の単元、すなわち非零定数になる。

divergence の式 `f_w=2l+y l_y, f_z=−2m−y m_y` を y=0 へ制限すれば
`l=αy,m=βy` となり、

    f=c+uy+vy²+3y(αw−βz), c≠0

を得る。α=β=0 なら Casimir 条件により F は y のみの多項式となり、
四次斉次化は二変数同次式なので複素数体上で分解する。
それ以外は `aα+bβ=1` として `z'=az+bw,w'=αw−βz` とすればよい。
行列式は1であり、α≠0 と β≠0 の二主開で全てを覆う。
この操作は Z を保ち、`l=y,m=0` とする。

Casimir 方程式は `F_z=0` と

    (c+uy+vy²+3yw)F_w−y²F_y=0

である。w の最高次数 j の係数 `A_j(y)` は `y²A'_j=3jyA_j` を満たす。
従って `A_j` は `y^(3j)` の定数倍で、全次数四以下から j≤1。
全解は `F=F0+kI`、

    I=y³w+cy²/2+uy³/3+vy⁴/4

となる。`I_y=yf,I_w=y³` を使い modular 恒等式の x 成分を計算すると
`aF_y+dF_w=ky(af+y²d)=kyF=−yF`。整域性と F≠0 から k=−1 が従う。
この段階で、同じ退化因子を持つという条件だけに証明を弱めていない。

`F=af+y²d` を z 微分し、`gcd(f,y)=1` と `deg a_z≤1` を使うと
`a_z=d_z=0` となる。y の零次・一次係数の比較は、正確に

    a=a0−a0uy/c−3a0yw/c+Ay²,
    d=(F−af)/y²

を与える。分子は y² で正確に割れ、商は二次以下である。
残る divergence から

    Pxz=xy+z(−y−a_y−d_w)+B₂(y,w)

であり、B₂ の二次以下六係数が全て自由に残る。
これで全十一変数の源が得られ、消した枝は既約四次と矛盾するものに限られる。

## producer・CAS を使わない算術再構成

実装は `code/classification_minimal_quartic_review.py`。
標準ライブラリの Fraction・疎多項式演算を用い、既存の独立 cofactor Ring を再利用した。
許す Laurent 逆元は明示的に c のみ。共通の全90係数 reviewer は変更していない。
有限 job `jobs/hyperplane_quartic_minimal_nilpotent_review.json` の gate を通過後に実行し、
0.0113秒、最大RSS 13,056 KiB で全検証が成功した。
結果は `reports/CLASSIFICATION_HYPERPLANE_QUARTIC_MINIMAL_NILPOTENT_REVIEW_20261107.json`。

PDE については、`y^i t^j` への作用を指数から独立に計算して全15列を生成し、
producer の行列の全係数へ照合した。
完全核候補 1,I を全行へ代入し、選択された13×13多項式 minor を
疎 Laplace 展開で実際に計算した。結果は非零定数14929920。
核の二つのベクトルは定数係数と `y³t` 係数で常に独立なので、全 c,u,v で rank13・核次元2となる。
標本点の rank だけからの推論ではなく、パラメータの例外値もない。

十一変数源では d を producer 表示から読むのではなく、上の Pf 恒等式から正確除算で再生成した。
全六括弧の一致、全四 strict Jacobi、全 divergence、実 Pf、
完全な modular 恒等式、全座標次数二以下を再展開した。

旧 MINIMAL_M の J+Q を `legacy/v07/minimal_thickness.tex` の六括弧から別に作り、
`(x,y,z,w)=(Y/6,X,−W/2,Z)` を作用させた。
`a0=−1/6,c=1/2,u=v=0,A=−1/2,B₂=−y²/2` における全六括弧と完全に一致した。
因子だけの一致や定数項を落とした一致ではない。

## 全既約源を一つの既知成分へ送る論理

`(C*)²×A⁹` は既約であり、この源の像閉包も既約である。
その像閉包を含む全 Poisson 空間の既約成分 C を一つ選ぶ。
像には上の実際の MINIMAL_M 点が含まれる。
`legacy/v07/minimal_thickness.tex` 定理minimal-fullの完成局所環は唯一の局所最小素を持つので、
この点を通る被約既約成分は MINIMAL_M 軌道閉包だけである。
従って C はその成分となり、源の全点が同じ閉成分へ入る。
GL・平行移動・非零 scalar・PGL の戻しは連結群作用であり、各成分を保つ。

全源の F は w の一次式で、係数 `−y³` と残る定数項は a0c≠0 により互いに素。
従って既約で、四次斉次化にも h 因子はない。
全パラメータで rank 四と単純 H が保持される。

ここで既存定理から使うのは「この点を通る被約成分の唯一性」だけである。
旧 nil 指数や形式局所環の再計算は実行していない。
独立に確認した源の完全性・全恒等式と、この既存定理の依存を区別する。
