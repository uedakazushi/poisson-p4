# 二H＋cubic三源のT38障害の独立確認と第三源の全T36障害

2026-11-08、PRIMARY E、別作者 `/root/log_manuscript`。
新operator、point、profile、minor、局所変形jetの延長又はCASは無し。
既印刷whole tensorの最低次数係数を手で比較する純証明である。
既知親への陽性包含、新成分の成分性又は全E完了は主張しない。
主稿/frontier/台帳/Git及び全producer bytesは不変。

## 1. T38の必要閉条件を全五parameterから独立に確認する

入力はtrue completed source
`reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json`。
座標(a,b,c,u,v)、p非零の全sourceは次の五つの非零entryを持つ。

\[
\begin{split}
 P_{ab}&=t_0(cv-2t_3u^2/p),\quad
 P_{ac}=t_1(bu+2t_2v^2/p),\\
 P_{bc}&=pbc+4t_2t_3uv/p,\\
 P_{bu}&=pbu/2+t_2v^2,\quad
 P_{cv}=-pcv/2+t_3u^2.
\end{split} \tag{1}
\]

全三correctionを含む。(p,t0,t1,t2,t3)の一般点でp*t2*t3が非零とし、
u,v非零、a自由に

\[
 b=-2t_2v^2/(pu),\qquad c=2t_3u^2/(pv) \tag{2}
\]

を代入する。Pbu/Pcv零、Pab/Pac零に加え、
bc=-4t2t3uv/p²だからPbcも零。
これは全parameterのliteral恒等式であり、t0/t1を零に置いていない。
(a,u,v)の三次元affine coneに属するraw zeroの稠密開を与えるので、
射影zero locusは曲面を含む。

普遍zero incidence
\(\{(\pi,q):\pi(q)=0\}\subset\mathcal P\times\mathbb P^4\)
は閉であり、基底への射はproper。
従って \(\dim Z_0(\pi)\ge2\) はPoisson parameter上の閉条件である。
既受理 `MPS_S61_CATALOGUE_20261103_v62.md` §2は、完成した全period/edge源のPGL飽和像閉包が
ラベルT38の全成分そのものであることを別途証明している。
よってこの閉条件はT38**全閉成分**に成り立つ。
特定のcoordinate cell又はperiod非零chartだけの障害ではない。

## 2. 三つの全lambda線のzero locus：producerの結論を別に手導出する

`TWO_LINES_PLANE_CUBIC_THREE_QUOTIENTS_SEPARATE_REVIEW.md` の全六磁気／actual inverseを入力にする。
以下では正規化後の全lambda源を扱い、正規化で反転した磁気係数の零補集合へ陰性を外挿しない。

第一・第二は

\[
 P=U\wedge\partial_z+V\wedge\partial_w+(-dw+\lambda h^2+hy)\partial_z\wedge\partial_w.
\]

第一の全anchorは
\(U=2x(x-h)\partial_x+y(3x-h)\partial_y\)、
\(V=-2xy\partial_x+(2x(h-x)-3y^2)\partial_y\)、d=x-h。
第二は
\(U=2x(x+h)\partial_x+3xy\partial_y\)、
\(V=2xy\partial_x+3(y^2-(h+x)^2)\partial_y\)、d=x。

h一でU=V零を直接解くと、第一は(x,y)=(0,0),(1,0)、
第二は(-1,0),(0,1),(0,-1)を尽くす。
第一の(0,0)、第二の(-1,0)ではd=-1なので、残る条件はw=-lambda、z自由で一本のzero line。
残点ではd零なので、第一はlambda零、第二はlambda=-1又は1の場合だけ全fiber planeがzeroとなる。

第三の全sourceは

\[
\begin{split}
 P_{xz}&=-hx,&P_{xw}&=3x(h-x),\\
 P_{yz}&=3(h-x)^2-3hy,&P_{yw}&=-9xy,\\
 P_{zw}&=3xz+3(\lambda h^2+y^2),&P_{xy}&=P_{h\bullet}=0.
\end{split} \tag{3}
\]

h一でPxz零からx零、Pyz零からy一、残りは3(lambda+1)。
従ってlambda=-1の時だけ、このchartに全zero plane x=0,y=1がある。
他lambdaにはh非零のzeroは無い。

h零かつx非零では、三つのraw四次元bivectorのPfaffianは非零定数倍のx四乗である。
標準pair順での係数は第一4、第二6、第三-9となるが、ここで使うのは非零性だけ。
raw rank四のPがEuler方向とのdecomposable外積になることはなく、projective zeroでない。
h=x零、y非零では、第一・第二のraw tensorはそれぞれ
-3y² partial_y∧partial_w、3y² partial_y∧partial_wである。
Euler wedgeの全条件はz零となり、それぞれ一本のlineだけを与える。
第三ではraw tensorは3y² partial_z∧partial_w、Eulerのy成分が非零なのでprojective zeroでない。
h=x=y零のvertex lineは全三源でzeroとして残る。

以上で

\[
\begin{array}{c|c}
 \text{sourceの開集合}&\dim Z_0(\pi)\\\hline
 \text{第一 }\lambda\ne0&1\\
 \text{第二 }\lambda\notin\{-1,1\}&1\\
 \text{第三 }\lambda\ne-1&1
\end{array}
\]

を全域の手比較から得る。
よって各開集合の全点はwhole T38の外。
元の全A6 sourceの飽和閉包の**親としてT38は三源とも不可能**である。
例外lambda、正規化磁気零又は別の既知親を除いたとはしない。
これは `TWO_LINES_PLANE_CUBIC_ZERO_LOCUS_T38_EXCLUSION.md` の幾何結論の独立確認である。
同producerに印刷されたPfaffian符号は上記標準pair順と反対だが、元証明も非零定数倍だけを使用しており、
zero-locusの場合分けと閉性には影響しない。旧bytesは変更しない。

## 3. Whole T36の正確なproper jet条件

`GA9_LINEAR_W_WHOLE_T36_QUADRATIC_PENCIL_EXCLUSION.md` §§1–3の閉incidence補題を入力にする。
四次元tangent空間Wのhomogeneous quadratic bivector Qについて

\[
 \alpha_Q=\iota_{E_W\wedge Q}\operatorname{vol}_W
\]

はEuler-horizontal cubic one-formである。
独立quadrics(F2,G2)のGrassmannianから
\([G2\,dF2-F2\,dG2]\) へのmapは全域でregular、properなので像R22は閉。
共通factor及びそのaffine coneの原点も含める。
「あるorder二以上のzero qがあり、alpha_Qがこの閉cone内」という条件も、
zero/first-jet incidenceとqのpropernessにより全tensor parameter上で閉じている。

true T36の指定vertexの全quadratic jetは、全period/edgeで
\(q_+=(p/5)Bu+t_2v^2\)、\(q_-=(p/5)Cv-t_3u^2\) により
\(\alpha_Q=(5/p)(q_-dq_+-q_+dq_-)\) と書ける。
従って**T36全成分**はこのproper marked-zero条件を満たす。
この節は既受理純補題を再利用し、新operator又はzero点探索をしていない。

## 4. 第三源の全候補qと、両vertex chartの全二次式

以下は第三(3)だけを扱う。§2により、lambda非-1のzeroは
\(L=\{h=x=y=0\}=\mathbb P(\operatorname{span}(z,w))\) だけ。
lambda=-1の追加zero planeのh一部分では、Pxz=-xのx微分が-1であり、
**全点でfirst jetが非零**である。
そのplaneの無限遠は同じLなので、全lambdaでorder二以上のzero候補は厳密にLだけである。

base b=(h,x,y)、EulerをEbとし、

\[
 U=-hx\partial_x+(3(h-x)^2-3hy)\partial_y,\qquad
 V=3x(h-x)\partial_x-9xy\partial_y.
\]

q=[z:w]=[1:kappa]で、projective chart z一、r=w-kappaを使う。
比のPoisson括弧を全raw tensorへ適用すると

\[
 P_{bb}=E_b\wedge U,\quad
 P_{br}=V-(\kappa+r)U-(3x+3\lambda h^2+3y^2)E_b.
\]

前者は全degree三、後者のdegree二部分だけを取るので

\[
 Q_q=(V-\kappa U-3xE_b)\wedge\partial_r. \tag{4}
\]

q無限遠[z:w]=[0:1]ではw一、r=zとし、同じ全比の式から

\[
 P_{bb}=E_b\wedge V,\quad
 P_{br}=U-rV+(3xr+3\lambda h^2+3y^2)E_b,
 \qquad Q_\infty=U\wedge\partial_r. \tag{5}
\]

いずれもzero/first jetは零で、磁気lambdaはdegree三以上に残っている。
一つの指定q又は有限kappaだけを扱っていない。
homogeneous q=[zeta:omega]について、radial項を除いた全式は非零scalarまで

\[
 W_q=\zeta V-\omega U\pmod{(\text{linear form})E_b}. \tag{6}
\]

従ってalpha_Qはbaseからpullbackした
\(\iota_{E_b\wedge(\zeta V-\omega U)}\operatorname{vol}_b\) の非零scalar倍。
Saito行列式det(Eb,U,V)は非零定数倍hxC3なので、どのqでもalpha_Qは非零である。

## 5. Invariant cubicはquadric pencilと両立しない

\(C_3=h^2y-(h-x)^3\) はC[h,x]上のprimitiveなy一次式なので絶対既約である。
全係数の手微分は

\[
 U(C_3)=-3hC_3,\qquad V(C_3)=-9xC_3
\]

を与える。従って全qのalpha_Qのfoliationは、base+rのP3におけるC3=0を保つ。
このcubic上でalpha_Qが恒等零なら、全cubic係数はC3の定数倍であり
\(\alpha_Q=C_3\sum c_i db_i\) と書ける。
Euler水平性から \(C_3\sum c_i b_i=0\)、従って全ci零となり、§4の非零性に反する。
よってC3一般点でfoliationは非零、C3のtangent planeに等しい。

仮にalpha_Qが独立quadric pencil(F2,G2)から来れば、
その有理第一積分J=F2/G2はC3上で定数である。
degree二の分母はirreducible cubic C3を含めないので、この制限は一般点で定義される。
char零より \(F_2-cG_2\) がC3で割れなければならず、degree二対三から零polynomialとなる。
これはF2/G2の独立性に反する。
したがって全q∈Lの非零alpha_Qは閉cone C22の外である。

§4は全order二以上zeroをLに尽くし、§5はその全qを除くので

**第三cusp接線・尖点割線の全有限lambda源はwhole closed T36の外である。**

lambda=-1にも追加zero planeの全simple確認を入れたため、T36についてはこの値を例外にしない。
T38については依然lambda=-1を除外していない。
元の正規化磁気I零へいずれの陰性も外挿しない。

## 6. 残るE義務と再開点

第三全A6の23次元閉像は、whole parentとしてT36/T38へは入らない。
T37へのactual whole-lambda曲線、又は別の全componenthood／全prior非所属がなお必要。
上記二つの陰性だけで新成分数を増やさない。
全既知親の他の陰性を使ってT37又は新23成分へ絞る場合も、それぞれの受理済みscopeを別に保持する。
新しい局所高次jet、固定frameのpoint列又はoperator列へ進む理由にはしない。
