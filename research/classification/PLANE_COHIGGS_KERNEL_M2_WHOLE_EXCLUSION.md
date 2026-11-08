# Co-Higgs rank二・M(p)零・m二：全sourceの既約quintic空性

2026-11-08、PRIMARY E、producerの純証明とONE有限Fraction証明書。
**producer全源証明完了、独立受理は別義務**。
新component、全co-Higgs、全Eの完成は主張しない。
1054主原稿／frontier／global／Git／旧producerは変更していない。

## 1. 命題と全入口

P4上の一般rank四Poisson構造がrank零plane Pを持ち、実退化五次Dが被約既約とする。
normal位数m_D=2、co-Higgs共通image O(1)、M:C3→sl2のrank二、
inclusionのpoint pでM(p)=0なら、この全sourceは空である。
全irregular/higher、全normalizerの逆、全rhoと零境界を含む。

全first symbolと全実逆は既受理
[kernel独立review](PLANE_RANK_ZERO_RANK_TWO_KERNEL_SEPARATE_REVIEW.md)
によるND、BSS(rho)、BN1、BN0に尽くされる。
全30のP1、全九のP2、Euler六商と全350元式を使う。
新しいlinear Jacobi operator、rank、minor、point又はjetを構成しない。

[全shear補題](PLANE_COHIGGS_SHEAR_EXACT_SOURCE_REDUCTION.md)により、
ND、BSS全rho≠−3（rho=±1,3の壁も含む）、BN1では
kerD1=actual shear六、D2単射なので、全finite inverse shear後にP1=P2=0。
元P0の実五次はそれぞれ

~~~
4bc(bz²−cw²),  −6bcw(2bz+cw),  2b²w(2bz+cw).
~~~

従ってこれら全源は空。残るBSS(−3)とBN0だけを以下で扱う。

## 2. 全七核から六shearと一つのNへの全逆

旧完全kernel Kは30×7、旧actual shear Sは30×6。
旧S証明書の同じ非零六minorのrow集合Iを使い、

~~~
A=(S_I)^(-1)K_I,   R=K−S A.
~~~

Rの一つの非零列を非零coordinate一で正規化してNを取る。
全30×7でR=N bをliteralに確認した。
N_I=0でN非零なのでNはshear外、K=S A+N bは完全な直和分解である。
全A、b、Nを結果reportへ保存した。新minor/rankは使わない。

P1=K alphaに対し、全六係数A alphaのfinite inverse shear

~~~
(a,b,c,z,w) ->
(a+s0*z+s3*w, b+s1*z+s4*w, c+s2*z+s5*w, z,w)
~~~

でP1=tNへ戻り、t=b alphaである。
逆は全sの符号反転、det一。全六sを独立変数として正逆両合成、
保存Sとの全weight一Euler商等式、全曲線の150係数push/inverse roundtripを確認した。
誘導される全P2の二次補正もfinite変換と一緒に動かしており、P2を先に零とは置かない。

### BSS(rho=−3)

保存順の七kernel座標をalpha0,…,alpha6とすればt=−alpha0/4。
正規化した唯一の残方向Nの非零上三角括弧は正確に

~~~
Nab=cw/6, Naw=8zw/3, Nbc=cz, Nbw=−zw, Ncw=−4z².
~~~

### BN0

t=−3alpha3/2。残方向は正確に

~~~
Nab=5(a−c)w/6,
Nac=(a+13c/6)z,
Naw=−2zw,
Nbc=8bz/3−5cw/6,
Nbz=−5zw/6,
Ncw=−2zw/3.
~~~

それ以外のN括弧は零。これはchosen pointのperturbationではなく、
両全七parameter sourceの六shear商を全て記述する一parameterである。

## 3. FULL九P2を旧D2 inverseで一意に復元

旧D2は350×9、全両源でrank九である。
同じ旧非零九minorのrow集合Jとcolumn順を使い、P1=tN上のweight二方程式から

~~~
P2=t²Q,  Q=−(D2_J)^(-1) Jac(N)_J
~~~

を得る。選択式だけを解いたことで余分なsourceを消してはいない。
任意の全Poisson solutionはこの九式を満たすので、D2単射性により必ずこのP2である。
残る元全350式を全て戻すのが次節である。

旧九方向の順序は(ab:z²,zw,w²; ac:z²,zw,w²; bc:z²,zw,w²)。
得たactual Qは

~~~
BSS(−3): Qab=−z²/18,
BN0:     Qab=59z²/108,  Qbc=5z²/6,
~~~

他の全Q成分零。
元の全九未知はこの解を導く前まで保持し、free P2を捨てていない。
旧核に食い違いがあれば停止するguardも実際に通過した。

## 4. 元の二行でt零が必要

元quadは両方(a,b,c,w)。全350 rowは旧独立checkerの
quad lex、degree四monomialの固定順である。

|全source|元row|空間monomial|全P0+tN+t²Q上の係数|
|---|---:|---|---|
|BSS(−3)|121|bz³|−32t²/3|
|BN0|122|bz²w|−16t²/9|

これはJacobiの射影Euler quadrupleの係数であり、strict三添字式との取り違えは無い。
全350、全parameter次数0..4をliteralに展開し、
次数0／一の全係数零、上記非零二乗係数を確認した。
したがって複素被約支持では必ずt=0。
形式nil指数や全scheme idealの計算は不要である。

t零でP2も零、元P0はBSSで−6bcw(2bz+cw)、BN0で2b²w(2bz+cw)。
いずれもproper factorを持つので、指定する既約五次のrank四源は空である。
actual PGL変換の逆を全higherへ戻しても、既約性とrankは不変。
全first-symbol逆を戻し、§1の他の全rhoとND/BN1も合わせて冒頭の全空性を得る。

## 5. 有限証明書と実際の実行

runnerは新しいdistinct caller
[classification_plane_cohiggs_kernel_pointwise.py](../../code/classification_plane_cohiggs_kernel_pointwise.py)。
旧別作者Fraction疎多項式の演算だけをimportし、producer／CASはimportしない。
旧shared runner又はcurrent契約は変更していない。

job [plane_cohiggs_kernel_pointwise](jobs/plane_cohiggs_kernel_pointwise.json) は
全入力hash、180秒／3GiB、全二source一回だけの成功・非閉鎖分岐を先に固定した。
最初のgateはtask_idが非activeのEで失敗し、計算は開始しなかった。
active E_reduced_nonnormalへ直した契約でgate PASS後、一回実行した。
実測0.136秒、最大RSS16,808KiB、全明示恒等式PASS。
結果 [CLASSIFICATION_PLANE_COHIGGS_KERNEL_POINTWISE_20261108.json](../../reports/CLASSIFICATION_PLANE_COHIGGS_KERNEL_POINTWISE_20261108.json)
は実行source全文、全N/Q、全K=S A+N b、全350残式、actual Euler Pfaffian全126と全有限逆を保存する。

|artifact|SHA256|
|---|---|
|runner|e631f9dc2814f58586125564aa1870840456310cead44f2a9bec4f4571128521|
|job|9206b3db4aea41908ed75f97e8cad5cd1cc8829c11edf1b6c927ab5a33cf1ff7|
|result|2fc17b6aef1025645c072c84faf67a4d1ccfe7caec07da370fa0e44b3e8158f3|

このproducer証明の独立受理、global claim採録と主稿統合はrootの別義務である。
rank M三tau零、rank M二M(p)非零の全八family及びhigh-m四source、
rank M一以下、O(0)、scalar、zero-symbol又は全Eへ空性を外挿しない。
