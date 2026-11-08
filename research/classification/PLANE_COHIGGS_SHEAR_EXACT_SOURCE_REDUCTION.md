# Rank零plane co-Higgs：全shear核からの点ごとの還元

2026-11-08、PRIMARY E、担当 log_manuscript。新しいM1、point、profile、
minor又はjetを作らず、既受理の全30/9 sourceと実shearを使う純証明である。
1054主原稿／E inventory／global／Gitは変更しない。

## 1. 全sourceの有限shear補題

rank零planeをz=w=0、baseを(a,b,c)とする。
既受理のcanonical Euler商で全tensorはP=P0+P1+P2、
P1は全weight一30方向、P2は全weight二9方向を持つ。
元全350式のweight一はD1(P1)=0、
weight二はD2(P2)+Jac(P1)=0である。

実unipotent linear transformation

~~~
G_C: x_i -> x_i+C_i0*z+C_i1*w,   z -> z, w -> w
~~~

の逆は全六係数で正確にG_{-C}、detは一。
作用はnormal weightを上げるので、P0のweight一変化は保存済みの30×6 shear S、
P1のweight一部は不変、変化はP2だけに入る。
Euler補正を含めたcanonical sourceを使い、誘導P2を捨てない。

**補題。** ker D1=im Sかつker D2=0なら、
全actual Poisson sourceはP0の実PGL軌道へ戻る。

実際、任意P1をS Cと書いてGの符号を選びP1を消す。
全finite transform後もPoissonだからweight二はD2(P2')=0、
したがってP2'も零。これはnormalizationのgeneric点又はformal jetの主張ではない。
全係数Cは任意、全global逆を持つ一回のlinear changeである。

## 2. M(p)=0、rank M二、m二の直接の適用

PLANE_RANK_ZERO_RANK_TWO_KERNEL_SEPARATE_REVIEW.mdの全actual normalizerを戻す。
ND、BSSの全rho（rho=-3だけを除く。rho=±1,3も含む）、BN1では、
全D1核六、実shear rank六、D2単射が既に全350行で独立確認されている。
従って全higher sourceが各P0へ戻る。
保存されたactual五次は

~~~
ND: 4bc(bz²−cw²),
BSS: −6bcw(2bz+cw),
BN1: 2b²w(2bz+cw).
~~~

全てproper productであり、指定する被約既約quintic源は全空である。
元rhoの全parameterは保存し、rho=-3を境界として黙って閉性へ入れない。
BN0とBSS(-3)はker D1七なので次の有限全源問題へ残す。

## 3. 残二whole sourceで必要な一回だけの有限判定

各残源の保存済みSの非零六minorと完全核Kを使う。
その同じ六rowで投影して

~~~
K = S A + N b
~~~

という全30行の恒等式を確認する。Nは残一次元のactual bivectorである。
新rank/minor又は別pointは要らない。
全finite shear後のP1はtNとなる。
既存D2の同じ非零九minorから一意な候補P2=t²Qを解き、
元全350式へP0+tN+t²Qをliteral代入する。

非零残式があれば各係数は単項式constant*t^j（j=2,3,4）なので、
複素被約支持はt零だけ。P0のproper因子で全空となる。
全残式零なら実whole polynomial曲線を得る。
その全actual五次にproper因子があれば全空、
因子が得られなければこの全一parameter源が次の正の認識義務であり、
計算成功をE完成へ勝手に読み替えない。

この全二源を一batch、180秒／3GiBで判定する。
旧D2核が単射という前提と食い違えば停止して全free P2を保存する。
新M1、追加anchor、GB又は高次jetへ進まない。

## 4. 他のco-Higgs源との射程

この補題は旧rank M三又はnonkernel源のshear-exact層にも適用できるが、
全例外壁／全D2 kernelを確認していない全familyを一括で空とはしない。
特にrank M三tau零の二源、M(p)非零high-m四源、
nonkernelの全二parameter壁は別の未処理義務として保持する。
今回の新結論は上記M(p)=0／m二のshear-exact層の全空性と、
残二源への完全な一parameter還元である。全E／分類主定理／新成分数を主張しない。
