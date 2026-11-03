# 全S61の42-normal基底上の有限5生成表示と残余の閉埋込み

2026-11-03。一般betaの既に指定した非空開で、全アフィン代数Aは**42個のnormal座標nの多項式環R上で、五元により有限生成**である。さらに既知五点を分離した残余B=A/(e)のnormal projectionは**閉埋込みSpec B→A42_n**となる。betaそのものの上での有限性、R-moduleの自由性、残余の空性、一般fiber長さ5は未証明である。

これは既存の42×42行列式E42を展開する計算や、次数1の全left inverse再試行ではない。線形normal因数化の**z*n cross部分だけ**に対する定数126×126逆行列が新しい入力である。

## 1. 126個全てのcross積を解く

凍結済み全134式は、普遍既知三平面の座標z3、normal n42で、五本F0..4と

    G_i=A_i(z,n)*n, i=0..128

に定数可逆行操作で書き換えられる。各G_iを

    G_i=sum_(k,j) C_(i,kj)*z_k*n_j + D_i(n)

と分解する。D_iは次数≤2で、定数項はない。**Cは129×126でQQ rank126**である。mod32003で選んだ126行の正方行列に対し、installed FLINTのfmpq演算で両方の逆行列積を全係数検査した。全逆行列係数の分母は32003単元である。さらに全126×945 normal係数で

    P=−C_selected^-1*D_selected,
    z_k*n_j=P_kj(n),  degree P_kj≤2                    (1)

をQQで検査した。他の三行は(1)を引くと三本のpure-normal二次式H(n)となる。全129行を失わない。

betaを変数に戻した場合も、保存された有限の126×126逆行列recipeを同じ普遍枠上で使える。中心のminorはp-unitなので、betaだけの開集合を局所化してp-integralな表示を得る。ここでは一般巨大係数をserializeしたとは主張しない。

## 2. 五生成の正確な有限性

(1)でmixed monomialを減らすと、全てのnを含む単項式は純n多項式になる。zの個数が減るため停止する。五本F0..4も

    f_l(z)+K_l(n)

に変わる。f_lは既知平面との真の交差の五本の二次式であり、K_lの定数項は零である。

f_lのz最高二次部分は六個の二次単項式中でrank5である。その一次倍15本の最高三次部分は三変数の全十個の三次単項式を張る。これら小行列のQQ逆行列と32003-unit分母を保存する。従って任意の純z三次単項式は、次数≤2の純z項とnを含む項へ減り、後者は(1)で純nへ減る。五本の二次rewriteにより二次単項式は一個だけ残る。

よって

    1, z0, z1, z2, 一個の二次単項式

が**実際の全AをR上で生成**する。これは自由性を仮定したborder-basis主張ではなく、元I内の正確なrewriteから得る有限生成の上界である。normal=0のfiberは既知E5そのものなので、少なくともその点では五元が必要である。

従来の17-base/1031生成表示とbaseが異なる。42-normal base上の五生成をbeta-only有限性と混同しない。p-unit逆行列により、この有限生成は凍結p-integral中心を含むbeta局所基底上でも成立する。

## 3. 残余のnormal projectionは閉埋込み

先行のglobal splitterよりA=E5×B、N=(n)A=(1−e)Aである。AがR有限なのでBもR有限である。origin m=(n)ではB/mB=0であり、有限性とNakayamaによりB_m=0。一方D(n_i)では(1)からz_k=P_ki(n)/n_iなので、R_(n_i)→B_(n_i)は全生成元を含み全射である。

従って有限R-module coker(R→B)の全局所化が零となり、**R→Bは全射**である。normal射影はスキームとして閉埋込みであり、その像はoriginと交わらない。これは実際に証明したR-有限性を用いたNakayamaであり、不適切なmodp空性からの推論ではない。

## 4. 全残余を保つnormal-only低次数criterion

三本のH、三組の互換式

    n_i*P_kj(n)−n_j*P_ki(n), k=0..2,i<j

および

    n_i²*(f_l(P_0i/n_i,P_1i/n_i,P_2i/n_i)+K_l(n))

を全i,lについて取ったidealをJ⊂Rとする。最後の表示は分母が消えた次数≤4の多項式である。生成元の数は3本の二次式、2583本の三次以下の式、210本の四次以下の式、計2796本である。

D(n_i)では互換式から全(1)が戻り、最後の五本から全F0..4が戻る。従ってJの各origin外chartは実際のBとスキーム同型である。Bは既に閉でoriginを避けるので

    I_B_normal = J : (n0,...,n41)^infinity             (2)

が正確に成立する。空性の具体的目標はJのある有限次数までのMacaulay spanがmの全冪を含むこと、すなわち(2)が単位idealとなることである。これはまだ計算していない。

## 5. determinant-free projective incidenceとの一致

元の45座標と[v]∈P41に対して、F0..4、A(q)v、n_i v_j−n_j v_iで定めたZはSpec Bとスキーム同型である。残余ではNが単位idealなので[v]=[n]が強制される。既知因子ではn=0だがAの42-column rankが完全なのでprojective fiberは空である。chart v_i=1でn=r*vとすると129式は(z3,r)に線形、残る五式はこの四変数に二次である。rは実際の残余上で自動的に単元となる。

このprojective incidenceはA45上でproperだが、beta上でproperとは主張しない。(2)はこの全chartを保持したscalar版である。

## 6. 再現と保存

`code/s61_affine_normal_finite5_20261103_v63.py` は元の全129線形normal行列からC,Dを再構築し、`code/s61_affine_normal_cross_inverse_20261103_v63.cpp` をinstalled FLINTでcompileして正確な検査を行う。入力と小plane逆行列、ログ、巨大出力のsha256を `reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/` に保存する。107MBの展開された126逆行列と126×945係数は数秒で再生成できるためignored workに置き、38MB圧縮物はrepositoryへ追加しない。

先のSymPy逆行列試行は90秒で打切り、完了を主張しなかった。native法は別の正確な有限線形代数である。全結果は元I、beta-onlyの開集合、あるいは明示したnormal saturationに関するものであり、全S61の一般五点尽くしはなお未解決である。
