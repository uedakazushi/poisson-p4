# rank零plane・normal冪零・b非零の全source非優越：独立査読

2026-11-07、主義務G。別作者による幾何と exact Fraction 算術の独立確認。
この査読の結論は **m_D=2、global-vector supersource、normal冪零、b≠0 の全像の次元が22以下**
という非優越だけである。既知成分へのmembership、b=0、全plane、全G又は全Eの完了は含まない。

## 1. 入力、有限job、算術の独立性

固定入力は次の三件である。

- `PLANE_RANK_ZERO_GLOBAL_VECTOR_NILPOTENT_FIRST_JETS.md`。
- `PLANE_RANK_ZERO_GLOBAL_VECTOR_LOW_PARAMETER_WEIGHT_BOUND.md`。
- `reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json`。

全hashと使用する独立算術libraryのhashは
`jobs/plane_global_nilpotent_nonzero_b_review.json` に固定した。
180秒・1536MiBのjobをgate PASS後に実行した。gateは数学を検証しない。
独立checkerは `code/classification_plane_global_nilpotent_nonzero_b_review.py`、
証明書は `reports/CLASSIFICATION_PLANE_GLOBAL_NILPOTENT_NONZERO_B_REVIEW_20261107.json`。
実測2.507秒、peak RSS 26,784KiBで全検査PASS。

checkerは先行**独立reviewer**のFraction(tau)、疎多項式、Schouten、Euler商の算術だけを再利用する。
producer及びCASをimportしない。producerの行列を入力としてrankだけを再計算したものでもない。
下記四つのT,Vからtensorと全operatorを別に構成して、保存された全entryと照合した。
可変runnerを参照せず、入力report内のimmutable producer bytesのSHA
`947b14874c9ebccc9f94164816958a57829ce3136800c7c5e6dfcd411ac5bf44`
を再検証した。元の二十一層packetの `all_passed=false` はそのまま保持している。
本review reportの `all_passed=true` は、表題のb非零全sourceだけを指す。

## 2. 全weight零条件とb非零のactual正規化

P=(z=w=0)、Euler_normal=z∂z+w∂w とする。
global-vector入口により、定数normal sl2行列が非零冪零である全weight零sourceは

\[
 \Pi_0=T\wedge E_n+V\wedge z\partial_w
              +(az^2+bzw+cw^2)\partial_z\wedge\partial_w
\]

である。T,Vは三つのbase座標のtrace-free linear fieldであり、Vのdivisorial zeroも許す。
base関数f,z,wのJacobiを直接展開すると
`z²([T,V]_field+bV)f+2czwVf`。
base二関数を含むJacobiは自動なので、V≠0ならc=0と `[T,V]_field=−bV` が必要十分である。
projective Euler wedgeで消えるscalar matrixはtrace零から零になる。
leading Pfaffianが非零であるための条件は
`det(Euler_base,T,V)≠0`、特にT,Vは独立である。

ここでb≠0とする。normal shear `w'=w+(a/b)z` はE_nとz∂wを保ち、
Kをbzw'へ送る。逆はminusのshearである。tensor scalar 1/bによりK=zwとなる。
全higher coefficientsも同じ可逆操作で移すので、first symbol以外を切り捨てない。

係数行列AのfieldをAxと定義するとfield bracketはmatrix commutatorの負であり、
actual base changeはAをgAg^{-1}へ送る。従って正規化後の行列条件は `[T,V]=V`。
`[T,V^j]=jV^j` のtraceをj=1,2,3について取るとtr(V^j)=0。
Newton恒等式からVは非零冪零であり、Jordan3又はJordan21だけで全sourceを尽くす。
b=0へこの論法を延長していない。

## 3. 全Jordan3、全Jordan21と二つの共鳴

J=E01+E12、G=diag(1,0,−1)とする。
V=Jに直すと `[T,J]=J` の全trace-free解は `T=G+xJ+yJ²`。
九未知数にcommutatorとtraceを課す線形写像のrankは7であり、二つの解parameterを残す。
actual centralizer change `g=exp(xJ+(y/2)J²)` はJを保ちTをGへ送る。
指数と逆は有限多項式である。独立checkerは全x,y上でg g^{-1}=I、Ad_g T=G、Ad_g J=Jを検証した。
よってJordan3はwhole orbit

`R: T=G, V=J, K=zw`

だけであり、leading determinantは `x2(2x0x2−x1²)` である。

Jordan21ではN=E01、B=E02、C=E21とする。
`[T,N]=N` の全trace-free解は

\[
 T=D_\rho+xN+yB+vC,\qquad
 D_\rho=\operatorname{diag}(\rho+1,\rho,-1-2\rho).
\]

同じ九未知数写像のrankは5、全解次元は4である。以下の恒等式を
rho,x,y,v,s,t,uの全多項式環で独立に検証した。

| actual change | xの像 | yの像 | vの像 |
|---|---|---|---|
| I+sB | x+sv | y−s(2+3rho) | v |
| I+tC | x−ty | y | v+t(1+3rho) |
| I+uN | x−u | y | v |

すべてNを保ち、逆はparameterの符号反転である。
rho≠−2/3,−1/3なら順にy,v,xを消す。
rho=−2/3ではv,xを消した後のy、rho=−1/3ではy,xを消した後のvを保持する。
diag(d,d,k)はNを保ちBを(d/k)B、Cを(k/d)Cへ送る。
全Laurent恒等式と逆を検証したため、残係数非零は1へ直り、零は別層に残る。

全解でleading determinantは
`x1²*(v*x1−(1+3rho)*x2)`。
rho=−1/3,v=0だけはm_D=2の外であり、m_D≥3の別義務として保持する。
非零vの共鳴sourceを捨てず、全m_D=2 sourceは次の三つとなる。

| model | T | V | parameter |
|---|---|---|---|
| NDIAG | D_rho | N | 全rho≠−1/3 |
| NB | D_{−2/3}+B | N | 一つのwhole orbit |
| NC | D_{−1/3}+C | N | 一つのwhole orbit |

それぞれK=zwであり、determinantは `−(1+3rho)x1²x2`、`x1²x2`、`x1³`。
以上は全centralizerによるactual正規化であり、chosen pointの計算から全源を推論していない。

## 4. 全30/9源、全350式と全minor wall

weight一のraw quadratic bivectorは36次元で、六つのEuler gaugeを除いた全30方向を保持した。
checkerはraw36に対して選んだ30方向と六gaugeを合わせたrank36を独立に確認した。
weight二はEuler gaugeを持たず、rawの全9方向である。
入力のsymbolic sourceの各parameterを別々に抽出し、この全基底と一致することも確認した。

四モデルの各々について全150 homogeneous coefficientから正しいEuler-contracted Pfaffianを計算し、
保存された実quintic及びleading determinantと一致した。
full projective Schoutenの全350係数を用いて
`D1=[Pi0,−]:30→350` と `D2=[Pi0,−]:9→350` を独立に構成した。
六つのactual linear shear（normal座標をbase座標へ足す変換）の微分もEuler商で全30×6へ構成した。
全entryの一致、各nonzero maximal minor、complete kernelの独立性、全行へのkernel代入零を確認した。

| model | rank D1 / nullity | rank D2 / nullity | actual shear rank |
|---|---:|---:|---:|
| R | 24 / 6 | 9 / 0 | 6 |
| NDIAG、generic rho | 24 / 6 | 9 / 0 | 6 |
| NB | 24 / 6 | 9 / 0 | 6 |
| NC | 24 / 6 | 9 / 0 | 6 |

固定三モデルR,NB,NCのD1/D2 minorはそれぞれ
`(1152,−1440)`、`(18,144)`、`(−8,−144)`。
NDIAGのgeneric minorsは

\[
\begin{aligned}
 m_1={}&-54\rho(\rho+1)(2\rho+1)(3\rho-1)(3\rho+1)^4(3\rho+2)^6,\\
 m_2={}&54\rho(\rho+1)(2\rho+1)(3\rho-2)(3\rho-1)
             (3\rho+1)(3\rho+2)^2(6\rho+5),\\
 m_s={}&-(3\rho+2)^2.
\end{aligned}
\]

この三minorの積と保存されたexception polynomialを独立に照合し、
その全一次因子の積へ定数を除いて完全分解したことを確認した。
従って未処理の非有理wallはなく、全wallは

`rho=2/3,−5/6,0,−1,−1/2,1/3,−1/3,−2/3`

だけである。それぞれで全350行を直接specializeし、別のnonzero minorとcomplete kernelでrankを再認証した。
全wallでD2はinjective、shear rankは6。
D1 nullityは−1/3と−2/3で7、それ以外で6である。
−1/3のdiagonal rowはm_D≥3なので今回の結論に数えないが、算術自体は検査済みである。
同じrhoの非零共鳴row NCは別に全sourceで検査済みである。

## 5. 全affine symbol次元とfull Poisson fiberの復元

normal非零冪零lineのchoiceは1次元である。
Jordan3のVのmatrix orbitは6次元、Tの全解spaceは2次元、元のKのa,bは2次元で、
全fixed-plane affine first-symbol baseは `B=1+6+2+2=11`。
Jordan21ではV orbit4、T解4なので `B=1+4+4+2=11`。
rho固定又は二つの共鳴whole orbitではbaseに少なくとも一つの条件が入りB≤10。
Vの非零scalarはV orbitに含まれ、別に二重加算しない。
これらは全first-symbol空間の次元であって、normal-form parameterだけを数えたものではない。

任意のfull higher solutionは `Pi=Pi0+Pi1+Pi2` と書ける。
最低weightのJacobiはD1 Pi1=0。
固定したPi1の上でPi2はinhomogeneousなD2方程式の解集合なので、空又はker D2のaffine torsorの部分集合である。
残る高weight Jacobiはfiberを縮めるだけである。
従ってfull higher fiberは `d≤nullity D1+nullity D2`。
shearを商にしてからdimensionを仮定する必要はなく、ここでは六shearを引いていない。

plane choice6とtensor projectivization−1を戻すと全projective ambient像は `B+d+5` 以下。
R及びNDIAGのgeneric部分は `11+6+5=22`。
NDIAGの全非除外wallはB≤10なので、rho=−2/3では `10+7+5=22`、他では21以下。
NB及びNCは `10+6+5=21`。
actual normalizationsは全higher fiberを可逆に運ぶから、以上は全b非零sourceを覆う。

よって **normal冪零・b≠0・m_D=2 の全global-vector sourceの像次元は22以下**。
未知G成分の一般元をこのbranchへ置くことは出来ない。

## 6. 正確な未処理範囲

b=0では `[T,V]=0` だけであり、Vが冪零とは限らない。
そのCartan、mixed、regular-nilcentralizer、common-kernel及び他normal型はこのreviewでは受理していない。
二十一層producer内の十六good層を一括で独立受理したものでもない。
後続one-cut packetによる二つの非線形cutもこのreportには含めない。
rho=−1/3の零共鳴rowを含むm_D≥3、その他のhigher-multiplicity source、全E帰属、全分類は未解決として保持する。
このreviewによって成分数の下界を変更しない。
