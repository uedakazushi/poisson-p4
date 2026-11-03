# 全normal base上のglobal quintic表示：独立検証

`S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63` を受理する。十分一般のbetaに対する全S61代数Aは、42個のnormal変数の多項式環R上で

\[
A\simeq R[h]/(G,r_0,\ldots,r_{133}),\qquad \deg_hG=5,\quad G\text{ monic},\quad \deg_hr_i\le4
\]

と表される。ここでglobalとは全normal affine空間上であり、未知qやcyclic determinantで局所化しない意味である。beta-onlyの既知の非零minorの開は使う。AのR自由性、beta上の有限性、一般五点尽くしは主張しない。

## 今回の独立係数検証

元134式の61変数データと保存T01中心点から、45個の中心化q変数による全134式を別の `Fraction` 疎多項式実装で直接再構成した。既存の `load_fiber` やproducerの多項式演算はimportしない。

保存された45×3のplane embeddingへ最初の五式F5を代入し、証明書の実際の五つのplane式f_iと全係数で一致することを確認した。三つのcofactor identityも別の疎多項式積で確認した。cofactorの次数は2、2、3、分母は全て32003-unitである。数え方を明確にすると、実際の非零係数項数は30、30、51である。producerの31、31、52という項数には、各関係で零である第五cofactorについて `Poly(0).terms()` が返す一つの零項が含まれる。この表示上の数え方は多項式恒等式や結果を変えない。

degree4/5のMacaulay行列を元のf_iから再構成した。選択した30×30、51×51 minorは独立な有限体Gauss消去で非零となり、QQ上の非零p-unit minorを与える。さらに、z1,z2のquartic表示p1,p2とmonic quintic gをcofactorのtargetから再構成し、全五つの (f_i(h,p_1(h),p_2(h))) がgで割り切れることとgのseparabilityを確認した。従って凍結中心のplane algebraが厳密に五次のprimitive algebraであることも確かめた。

全134の代入と除法には、producerのunivariate係数辞書と異なる、**46個の記号を平らにした疎多項式**を用いた。n0,...,n41、S1、S2、Hの45記号とhを区別し、全q代入を元134式へ行った。保存された全134組について

\[
F_i(q(h,n,S_1,S_2))=q_i(h,n,S_1,S_2,H)G+r_i(h,n,S_1,S_2,H)
\]

を全QQ係数で独立に確認した。元代入のh次数は8以下、全remainderは4以下である。これはplaceholder環の恒等式であり、標本代入や消失補間ではない。producerの新しい出力先での一回の再実行も、全証拠とbyte一致した。

## 一般betaへの有限なrecipe

generic体K上の既知五点のaffine spanはK-rationalな3-planeであり、plane frameとoriginはK上で選べる。originをgenericなT01 rootと仮定しない。三つのfree座標を固定してfrozen中心に一致するframeを選び、分母の非零beta開に制限する。

既受理のleading quadratic/cubic spanは、このplane algebraの長さを5以下に保つ。一方、実際の既知五点を全て含むのでgeneric長さは5である。hがprimitiveである条件は凍結中心で非零な5×5 determinantであり、同じbeta-only開で (1,h,\ldots,h^4) は既知plane algebraの基底になる。従ってp1,p2とgはrootを選ばずK上で定まる。

degree4以下の多項式空間は35次元で、plane algebraへの評価のkernelは30次元である。選択したdegree4 Macaulay minorのrank30が保たれる開では、そのimageはkernel全体に等しい。同様にdegree5では56−5=51次元となる。このため関係z1−p1、z2−p2、gを、それぞれ次数2、2、3以下のcofactorで表す同じ有限逆行列recipeが一般betaでも成立する。有限体の消失からQQへliftしているのではない。

## normal correctionsは独立parameterではない

既受理の全cross identity (z_kn_j=P_{kj}(n))、\(\deg_nP\le2\) を使い、normal因子のある単項式からzを一つずつ消す固定順のRedを取る。Pに定数項がないのでnormal因子は消えず、zの個数が減るため有限回で停止する。元五式は (f_i(z)+K_i(n)) に変わり、\(K_i(0)=0\)、\(\deg K_i\le2\) である。

plane cofactorを全Aへ持ち上げると、厳密に

\[
z_1=p_1(h)+S_1(n),\quad z_2=p_2(h)+S_2(n),\quad g(h)=H(n)
\]

となる。各correctionは \(-\operatorname{Red}\sum B_iK_i\) という固定されたR内の多項式である。cofactor内のzを全て消すとnormal次数は最大で元の総次数までなので、次数上界は4、4、5となる。このrecipeを巨大な42変数展開として保存したとは解釈しない。

placeholder環の全134除法恒等式へ、これらの実際のnormal多項式を代入できる。h=z0、n↔n、z1,z2を上式に送る写像により、提示したquotientとAの間に二つの写像がある。元134式の代入恒等式が一方のwell-defined性を、A内の正確なcofactor identityが他方の合成が恒等であることを保証する。従ってsupportやnilpotentを落とさないglobal algebra同型である。placeholderを自由な別source parameterとみなさないことが、この段階で必要である。

## 5×670の実module表示とdegree85の回復

\(C=R[h]/G\) はmonic性によりR自由rank5である。全AはCの商であり、そのkernel idealは、R-moduleとして \(h^jr_i\)、(0\le j\le4\)、(0\le i\le133\) の670列で生成される。従って \(R^{670}\to R^5\to A\to0\) は正確な有限表示である。複数の作用行列の非可換な閉包を仮定する議論ではなく、全ては単一のcompanion matrixの多項式作用である。AそのものがR自由とは限らない。

既受理のnormal commutator packetでは、実際のknown-five idempotent eはnormal多項式D84で、\(\deg D84\le84\)、\(D84(0)=1\) である。従って (1-D84=\sum c_i(n)n_i)、\(\deg c_i\le83\) と固定順で因数を一つ抜ける。cross identityにより

\[
z_k(1-e)=\sum_i c_i(n)P_{ki}(n)
\]

がA内で成立する。残余Bでは1−eが単位元なので、z_kは次数85以下のnormal多項式として回復する。この有限recipeは元のidempotent分解とcross identityから従い、展開した巨大determinantを新たに計算したものではない。Bの空性も意味しない。

独立codeは `code/s61_affine_global_quintic_review_20261103_v63.py`、receiptは `reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63_REVIEW.json` である。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/s61_affine_global_quintic_review_20261103_v63.py \
  --replay work/global-quintic-producer-replay \
  --output work/global-quintic-independent-review.json
```
