# generic B18の全actual N完成環：Nakayamaと81関係の独立検証

凍結済み `AMBIENT_B18_N_FULL_COMPLETION_20261103_v63` を受理する。ここでNは被約なactual target成分であり、source ramification scheme、部分射影像、actual CA/CB交叉Dとは区別する。対象はgeneric B18の開軌道の代表点における完成環である。

## 有限jetから全idealへの正確な移行

既受理のactual imageは \(A=R\oplus It\subset B=R[t]/(t^2-\nu)\) で、Rはactual target座標から作ったregular28次元底である。別に受理済みの全critical quotientとminimal Cohen presentationから、\(J=I+(\nu)\) の **全formal ideal pair** は、νを固定して

\[
R=\mathbf C[[g_{19},\nu,s_1,\ldots,s_6,w,v]],\quad
J=(\nu,s_1,\ldots,s_6,wv,v^2)
\]

となる。従って \(\dim_\mathbf C J/\mathfrak m_RJ=9\) である。この入力は、単にsource critical idealとIの生成数を同一視する仮定ではない。

独立受理済みのactual odd2jetは、I内の七つの元がJの全七一次classを張ることを示す。νの一次項もそのspan内にある。full R firstjetが零となるように七元を引いた別の二元は、Jの一次部分のkernelに含まれる実sourceの二次元plane上で \(-4wv,-6v^2\) という独立な二次制限を持つ。

任意の \(\mathfrak m_RJ\) の二次部分は、linear coefficientとJの一次部分の積の和なので、このplane上で零となる。従ってこの二元は **J/\(\mathfrak m_RJ\)** の残る二classを張る。これはI/\(\mathfrak m_RI\)だけの生成数下界から等号を推測する議論ではない。

よって \(J=I+\mathfrak m_RJ\) であり、有限R-module J/IへNakayamaを適用して

\[
I=J,　\nu\in I
\]

が全次数で従う。modulo νの座標変換をRへ持ち上げてν自身を保てるので、以後のrank2 relationも \(t^2=\nu\) のままである。Torや未証明のflatnessを使わない。

## 全actual完成環と全kernel

従ってnormalizerとactual imageは

\[
B=\mathbf C[[g_{19},w,v,s_1,\ldots,s_6,t]],
\]

\[
A=\mathbf C[[g_{19},w,v,s_1,\ldots,s_6,
 t^2,s_1t,\ldots,s_6t,t^3,wvt,v^2t]]\subset B
\]

となる。

\(b=(\nu,s_1,\ldots,s_6,wv,v^2)\)、\(Y_i=b_it\) とすれば、全kernelは36本のTaylor pair syzygyと45本の \(Y_iY_j-\nu b_ib_j\) から生成される。全81式を独立に作り直して保存式と比較した。さらにproducerのeliminationとは別に、18変数のsource ringからregular polynomial ring \(\mathbf Q[w,v,s_1,\ldots,s_6,t]\) への写像のkernelをMacaulay2で直接求め、同じ81式のidealに等しいことを確認した。coreの次元も9である。

完全性の代数的理由も確認した。productsで任意の式を \(R+\sum RY_i\) に戻すと、基底1,tの独立性から残るkernelはmonomial ideal Iのfirst syzygyだけとなる。その全first syzygyはTaylor pairで生成される。finite R-moduleの完成により、同じ記述がpower series完成環にも成立する。

## 数値、不良 locus、normalization

\(A=R\oplus It\) と \(R/I=\mathbf C[[g_{19},w,v]]/(wv,v^2)\) により、coneのdim、edim、depthは28、37、20である。非零tensorのscalar因子を除くprojective版はgが18個となり、27、36、19である。

\(\mathfrak m_A=\mathfrak m_R\oplus It\) のodd productは \(\nu I^2\subset\mathfrak m_R^3\) に入るので、実際に

\[
\mathfrak m_A^2=(\mathfrak m_RA)\mathfrak m_A
\]

となる。従って \(\mathfrak m_RA\) はreductionであり、finite R-rankが2であることからHilbert--Samuel multiplicityは2である。AをR自由とは仮定しない。

Aはregular domain Bのsubringで、tは \(t^3/t^2\) としてfraction fieldに入る。従って完成環自体がdomainであり、embedded associated primeはない。conductorはB内でIB、A内で \(I\oplus It\) である。\(B/A=(R/I)t\) が持つembedded nil classを、A自身のembedded primeと取り違えない。

ν=t²とνt=t³がAに入る一方、tは入らない。\(B=A[t]\) はelementary subintegral extensionであり、finite universal homeomorphismである。Bはnormalなのでseminormalizationでもある。したがって、このformal germの中では、任意の接触次数の二つの相異なるgeometric normalizer点やreduced curve arcが同じactual像を持つことはない。非被約fiber長2はなお存在する。

conductorの被約台は \(\mathbf C[[g_{19},w]]\) というsmooth20次元coneである。その外ではA=Bがregular。conductor上の局所化ではIのheightは少なくとも2であり、principalにはならない。もしAがそこでCMなら、regular R上のfinite moduleとしてfreeとなり、直和因子Iもrank1 freeとなって矛盾する。従ってこのformal germ内のSing、Nonnormal、NonCMは全てconductor supportである。これをglobal N全体の不良locusの分類とは呼ばない。

この結果はfullNを同定するが、CA/CB個々の完成型やDの完成型を決めない。先に予定されたpaired-flag三次計算を依存関係に加えず、実actual odd2jetと全critical quotientから直接証明している。

独立codeは `code/ambient_B18_N_completion_review_20261103_v63.py`、receiptは `reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json` である。

```sh
. /workspace/cloud-setup/activate.sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_B18_N_completion_review_20261103_v63.py \
  --work work/B18-N-completion-independent-replay \
  --output work/B18-N-completion-independent-review.json
```
