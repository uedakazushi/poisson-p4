# NEW source 座標を戻す完全な有限有理演算の手順

2026-11-04。受理済みのshape次数10とZNのprimitivityから、sourceの五座標へ戻す
**具体的な有限operator recipe**を固定する。巨大なgeneric adjugateの係数展開や
新しいInv14のsource展開は不要である。唯一の数値検査は、既存Krylov証明書の
同じ100成分の行列に対する二つの線形方程式で、新しいβ点を評価していない。
actual target体との等式の証明は別のfield-recovery companionへ渡す。

## 入力とliteral polynomial data

target側で使う関数を

\[
R=(p_3+4p_4)/p_0,\quad Z_i=p_i/F_4^2\ (i=0,1,2),\quad
\eta=N_{12}/F_4^3
\]

とし、底をK=QQ(R,Z0,Z1,Z2)とする。primitive GはK上の既約10次である。
`recipe.json.gz`には次を参照だけでなく全係数で埋め込んだ。

* Gの891項。列はinteger,s,r,R,Z0,Z1,Z2。
* linear PRSのC1=A(r)s+B(r)の全constant 10,314項とlinear 10,698項。
  係数の列はinteger,Z2,Z1,Z0,R,rであり、**constantが先**である。
* 保存済みのN79,N41,Fn,Jの全有理係数。
* 元のquadratic C,Dの全32×136係数、136 monomialの順序、Omegaと厳密な逆行列。
* 新しい有理な軌道代表の全16βと、従来の1032-node source circuitの全入力・出力。

各元ファイルと受理済みgateのSHA256はreportに固定した。約839KBの非圧縮JSONであり、
未知の係数を後で探す命令は含まない。保存したprogramは、以下のfield演算・多項式剰余・
行列式だけを使う。

## Mの有限構成とrの回収

L=K[r]/(monic G)、基底を1,r,…,r9とする。C1のA,BをLへ写す。
s次数d以下のPに対し

\[
H_d(P)=A^dP(-B/A)
       =\sum p_{ijk}r^iR^kA^{d-j}(-B)^j
\]

と定義する。右辺の有限和を使うためAの除算は要らない。d0=r(r−1)²として

\[
U=Z_0d_0H_1(J)H_3(N79)+25Z_0^2H_2(Fn)H_2(N41),\qquad
V=20d_0H_2(Fn)H_1(J)^2.
\]

すべてmonic GによるEuclidean remainderで基底に戻す。
既受理のbalanced identityはη=U/Vを与える。Mの第j列をUʲV⁹⁻ʲの係数vector
(j=0,…,9)とする。Δ=det Mは受理済みの同一old point/32003の証明書で非零である。

\[
a=\operatorname{adj}(M)e_r/\Delta,\qquad
b=\operatorname{adj}(M)e_1/\Delta.
\]

ここでe_r=(0,1,0,…,0)、e_1=(1,0,…,0)。adj(M)の各成分は
(-1)^{i+j} times「row j、column iを削除した9×9行列式」であり、行列式は
通常のsigned permutation sumとする。従ってこの記法は完全に指定された有限有理演算で、
generic linear systemの未知の存在だけを主張するものではない。

両等式をLで読めば

\[
r=V^9\sum a_j\eta^j,\quad1=V^9\sum b_j\eta^j,
\qquad
\boxed{r=\frac{\sum a_j\eta^j}{\sum b_j\eta^j}}.
\]

generic LではV≠0であり、分母はV⁻⁹だから非零である。special algebraでVの逆を
仮定する議論ではない。

数値検査はKrylov原証明書の同じ保存Mだけを読み、Fp上で二つのsolveを行う。
`old_matrix_solves.json`にa,bの全20成分と、M a=e_r、M b=e_1の全20係数一致を
保存した。これは新しい点・新しいprime・新しいtarget ratio評価を含まない。

## s,q,Bの回収と開条件

\[
s=-B(r)/A(r).
\]

ここでのB(r)はlinear PRSのconstant coefficientで、source座標Bとは区別する。
A(r)≠0は既受理の3/5有理mapのbirational性とlinear C1により保証される。
次を順に計算する。

\[
T=(r-1)R-2r,\quad L_q=8Rr(r-1)-6r^2-11r+1,
\quad N_q=(s-1)L_q-(r-1),\quad q=N_q/(40rT).
\]
\[
h=s^2-25rq^2,\quad F=100r(r-1)q^2-4rs+s^2+4r,
\]
\[
D=2359296000r^4(r-1)q^2(-s+5q+1),\qquad
\boxed{B=Z_0hF^2/D}.
\]

これらは既存の両方向Möbius回収とcompressionの式をそのまま使う。sourceの既知dense
chartに加え、F4,p0,Δ,A(r),V,h14等を可逆にする。境界全域へregularに延長した、
properな逆写像を得た、とは主張しない。

## h14を展開しない有限演算とxi

普遍不変式

\[
J_{14}=\{\operatorname{tr}(C^4),\operatorname{tr}(C^3D)\}
=\nabla_\beta T_0^{\mathsf T}\Omega^{-1}\nabla_\beta T_1
\]

を固定する。符号は元のΩによるこの順であり、Inv14 basisの第一列と一致する。
その普遍非零性は受理済みminor7584から分かり、NEW定義次数下界26によりNEW上でも
消えない。sourceへ代入してからshape変数を微分する方法は使わない。
全16 ambient β方向の微分を、埋め込んだquadratic係数から先に作る。

\[
dT_0=4\operatorname{tr}(dC\,C^3),\qquad
dT_1=\operatorname{tr}\bigl(dC(C^2D+CDC+DC^2)\bigr)
        +\operatorname{tr}(C^3dD).
\]

producer内の`j14_operator`はこの演算を任意の係数体上で実行できる有限関数である。
今回はsymbolic sourceを与えて実行・展開していない。

有理軌道代表packetの全16式βbar(xi,r,s,q,B)を埋め込み、

\[
\boxed{h_{14}(r,s,q,B)=J_{14}(\bar\beta(1,r,s,q,B))}
\]

と指定する。根の選択や新しい消去計算は要らない。新packetを経由しない代替も保持した。
既存1032-node circuitへu²=1600r、X=s/(48u)、c=Bu/40、W=(q+1)/3072を代入し、
quadratic étale algebraでJ14を計算してtrace/2を取る。既知tauがSp4作用で実現されるため
J14はtau固定で、同じshapeの有理関数h14になる。これは全係数の巨大展開を要求しない
有限の四則演算の別表示である。

全16source係数で検証済みのcocharacter scalingにより、

\[
F_4=\xi^2f_4,\quad f_4=(2/405)F,\qquad J_{14}=\xi^7h_{14}.
\]

h14≠0なので、最後は

\[
\boxed{\xi=(J_{14,\mathrm{target}}/F_{4,\mathrm{target}}^3)
                    f_4^3/h_{14}}.
\]

これで全五source座標(xi,r,s,q,B)の有限有理recipeが揃う。full actual invariant-fieldの
指数1への帰結は別のproof companionに整理する。NEWの定義方程式、全orbitの分類、
境界上での一意性、full hypersurface自体の有理性をこのrecipeから主張しない。

## 保存と再現

* `code/q_new_family_constructive_inverse_recipe_20261104_v63.py`
* `reports/Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63/recipe.json.gz`
* 同じ場所の`old_matrix_solves.json`、`report.json`

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_constructive_inverse_recipe_20261104_v63.py \
  --output work/constructive-inverse-replay
```

新規output先を指定する。上限60秒/2GB/32MiB、実測0.358秒。
generic symbolic係数展開、new point、norm、source26は実行していない。
