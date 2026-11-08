# cusp-cone CI22 の全二つの五源：projective Casimirと全rank四空性

2026-11-08、PRIMARY E-FINITE-SURFACE。
**全cone_cusp tangent/conormal Poisson sourceはrank二以下。指定rank四E枝は全空。**
Gの≤23上界をmembershipと読み替えたのでなく、全actualD零とwhole Casimirを新しく確認した。

## 1. 全支持の被覆と零補集合

surfaceは

\[
 S=V(q_1,q_2),\quad q_1=2x_0x_2+x_1^2+x_3^2,\quad q_2=2x_1x_2+x_3^2,
\]

vertexはx4方向である。
`NONNORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md` §7の全full16-sourceを用いる。
元式p13*p15=0、両方零の全源ではactualDが恒等零。
残る全二openは、それぞれ次の完全な五parameter線形sourceへ入る。

| case | 自由係数 | 他の非零座標 |
|---|---|---|
| I、p13=h≠0 | p6=f,p9=b,p11=g,p12=a,p13=h | p4=2a,p8=−f,p14=−b |
| II、p15=h≠0 | p3=f,p10=b,p12=a,p14=g,p15=h | p1=2a/3,p4=4a/3,p7=(−12b+9f)/2 |

表に無いp座標は零。全inverse caseと全350式の整合性は先行別作者receiptに含まれる。
p13=p15=0をgeneral kernelの特殊化で代表させず、元の全zero補集合として残す。
以下はa,b,f,gの零／非零を制限しない。

## 2. case IIはraw decomposable

保存された全150 sourceをmap IIへ戻すと、Pi_ij=0 (i,j≤3)。
従ってΠ=V∧∂4であり、射影化してもrank≤2。
Vの全四係数は五parameterに線形で、reportに全entryを保存した。
hidden higher source又はextra conormal sectionを落としたのではなく、full16 kernelからの全mapである。

## 3. case Iのwhole projective Casimir

Π=hP_h+(aV_a+bV_b+fV_f+gV_g)∧∂4 と書ける。
P_hの非零上三角entryは

\[
\begin{aligned}
 (02)&=x_3(x_1-2x_2), &(03)&=-(x_1^2+x_3^2),\\
 (12)&=-2x_2x_3, &(13)&=-x_3^2, &(23)&=-2x_2^2,\\
 (04)&=-x_3x_4, &(14)&=-x_3x_4, &(34)&=x_2x_4.
\end{aligned}
\]

V_iはx0,…,x3方向の全四成分を並べると

\[
\begin{aligned}
 V_a&=(x_1^2,0,2x_2^2,x_2x_3),\\
 V_b&=(-x_1x_3,-x_1x_3,-2x_2x_3,x_1x_2-x_3^2),\\
 V_f&=(0,(x_1-x_0)x_3,(x_1-x_2)x_3,x_0x_2-x_1^2),\\
 V_g&=((x_1-x_2)x_3,-x_2x_3,0,x_2^2).
\end{aligned}
\]

直接微分すると、両q_jについて

\[
 V_a(q_j)=2x_2q_j,\quad V_b(q_j)=-2x_3q_j,\quad
 V_f(q_j)=V_g(q_j)=0.
\]

従って四つのvertical bivectorはq1/q2をstrictにCasimirとする。
P_hについては、matrix-vector convention P#dq=P·grad(q) で

\[
 q_2P_h^\sharp dq_1-q_1P_h^\sharp dq_2
                    =2x_1x_3q_2 E.                  \tag{1}
\]

全五成分の多項式恒等式である。
ゆえに全h,a,b,f,gについて、右辺がEuler方向で消えるため、
**非定数の射影有理函数q1/q2がΠのCasimir**となる。
rank四の非退化開上ではbivectorが可逆なので、非定数Casimirは存在出来ない。
よってcase Iの全源もrank≤2である。
strict Poisson lift又はnormalization-flatを仮定していない。
raw contractionがEuler倍であることがprojectiveの正確な条件である。

## 4. actual Pfaffianの全係数readback

有限E契約をgate後、既存の全350／126配列を表の全二mapへ代入した。
さらに全150sourceからsigned Euler Pfaffianを直接構成して比較した。
両caseで全350 Jacobiが恒等零、全126 actualD係数も恒等零。
**0.161260秒、RSS51,460KiB、全exact PASS**。
case IIのdecomposable性、case Iの(1)と同じ結論を、元actual五次でも確認した。
source parameterの選択点や係数の一般値からの外挿ではない。

両openと既済のp13=p15=0補集合を合わせると、全complex Poisson閉点はrank二以下。
従ってcone_cuspを非正規曲面支持に持つrank四・実既約被約quintic枝は全空。
新componentの増分は零。rank二以下は既存の全支持義務へ戻し、
他のCI22 model、nil32 rank二、singularM1又は全Eへ結論を広げない。

実行契約は `jobs/ci22_cone_cusp_ranktwo_whole_factor.json`、
結果は `CLASSIFICATION_CI22_CONE_CUSP_RANKTWO_WHOLE_FACTOR_20261108.json`。
全map、全十括弧、全残差、inline source、入力SHAを保存した。
旧code、旧report、currentjob、主稿、global台帳、Gitは変更していない。
