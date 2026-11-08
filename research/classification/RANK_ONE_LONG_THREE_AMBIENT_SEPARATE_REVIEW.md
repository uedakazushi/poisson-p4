# rank一long三候補：三つの同一点全ambient独立検証

2026-11-08、quartic_review、PRIMARY E_multiple_hyperplane。**全三件PASS**。
producerのSymPy/general ambient engineをimportせず、既独立の標準Fraction疎多項式／行列演算で再構成した。
各々一つの事前宣言点だけを使い、新point/profile/minor探索をしていない。

## 元whole sourceから点を作り直す

座標(h,x,y,z,w)、r=y−x²/2、整数scale六。
入力g,U,H,Fからd=divUを微分で求め、U(g)=dg、(U+1)H=Fgを直接確認した。
非零定数g二例ではFを(U+1)Hとして独立に作り直し、
g零SCでは全六係数の自由quadratic Fをそのまま保持した。

|source|g|U|H三次部|
|---|---|---|---|
|CONST_G_MU1|1|−z∂z+w∂w|z²w|
|CONST_G_MU13|1|−z∂z/3+w∂w/3|z³|
|G0_SCALAR_CUBIC|0|−(z∂z+w∂w)/3|zw(z+w)|

定数g二源のH二次以下とSCのFは、全く同じ宣言値
1+2z+3w+5z²+7zw+11w²。
元source(19)の全六括弧を計算してdegree二斉次化し、六倍した。
保存150vectorを点の定義として使っていない。
得た全150係数はanchor／reportと一致した。
全十strict Jacobi、全350projective Jacobi、modular(0,6h,6x,0,0)、
全五Pf、全126signed quintic係数も全て一致した。
actualDはそれぞれ元sourceの−36*homogenization_degree5(rg+H)である。

## 全行列、元の選択minorと完全核

元全150方向Bに対しJ(P+B)−J(P)−J(B)を用い、各全350×150微分を独立に作った。
全Euler25列、全GL25列とscalarも元Pから作り直した。
全entryをproducerの保存matrixへ照合し、全て一致した。
保存した同じminorの行列式をFractionで再計算し、
全complete kernelの元全行消滅とfree-coordinate単位行列を検査した。

|source|rank D|rank Euler|rank actual51|その点のPGL軌道|
|---|---:|---:|---:|---:|
|CONST_G_MU1|100|24|48|23|
|CONST_G_MU13|100|24|48|23|
|G0_SCALAR_CUBIC|100|24|46|21|

全350行でD*actual51=0も確認した。
新minorを探したのではなく、既選択の下界と完全核による上界が一致する。

## 三つの全成分性に必要な最弱の矢印

各whole sourceのPGL飽和閉包Zが既約・正確に25次元、
全prior130及び相互に異なるという全源純証明は別作者の凍結入口reviewを固定入力とする。
今回の点はいずれもその元全sourceに属する。

全150 lift係数空間で、Zのaffine coneのEuler24逆像は50次元。
各点の全微分rank100によりZariski tangentも50次元なので、
局所次元と接次元が一致し、局所環はregular、全空間成分は一意となる。
projective Euler商後は25次元で、閉既約Zはその全成分に等しい。
従って三つとも新しい25次元全Poisson成分である。
generic sourceの下界を選択点のorbitだけで代用していない。
各点のorbit値は上表の固定点計算で、次節の全source上界なしにgenericへ流用しない。

## Generic orbit補足の独立純読解

定数g二whole源は全域で非零modular場Xを持ち、L_XP=0。
従ってprojective stabilizerのLie次元は一以上、全source orbitは23以下。
同点の既非零M48minorの開性でgeneric orbitは正確に23、商は25−23=2次元。

SC全源ではq0を中心として、Xに加え
(z−z0)∂y、(w−w0)∂yが全F／全binary cubic HでPoisson場となる。
一般にell∈span(z−z0,w−w0)についてxy成分のLie derivativeは
−ell/3−Uell=0。他の成分はy非依存とPzw=0により零。
三場はprojective linear Lie場として独立で、Eulerにはならない。
従って全source orbitは21以下、同点M46によりgeneric orbitは21、商は四次元。
これらの上界はPGL飽和と閉rank条件でも保たれ、generic開下界と両立する。
全特別点のorbit次元や全軌道分類は主張しない。

## 有限実行と範囲

fresh job `rank_one_long_three_ambient_separate_review.json` を実行前gate PASS。
予算180秒／2GiB、実測**2.035049秒、最大RSS21,816KiB、全三件ALL PASS**。
全input hash、inline program全文、既minorの正確なdetと核検査結果を
`CLASSIFICATION_RANK_ONE_LONG_THREE_AMBIENT_SEPARATE_REVIEW_20261108.json` に保存した。
producer／旧source／旧証明書／台帳／Gitは変更していない。

whole十七source routingの他の未認識族をここだけで閉じず、
これら三sourceの全境界包含と新成分採録は全源証明・統合側の別責任と区別する。
全E、総成分数上界又は全Poisson分類の完了は未主張。
