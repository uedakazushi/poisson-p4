# 残三つのreduced1+1+3源：同じ旧plane入口への正確な戻し

2026-11-08、PRIMARY E、log_manuscript。純手導出、独立確認前。
新operator、normal point、profile、minor、jetは作らない。既accepted irreducible-D proofの適用仮定を明示するための新入口であり、全既知親への陽性包含ではない。

## 1. 二simple Hの交差にあるrank二plane

`TWO_LINES_PLANE_CUBIC_THREE_MAGNETIC_QUOTIENTS.md` の正規化後の全I,J二源を使う。
planeはS=(h=x=0)、base(a,b,c)=(y,z,w)、normal(u,v)=(h,x)。
node-flexとcusp-flexではSに沿うD位数は三、先頭は非零定数×a²u²v。
S上の実Poisson zero cubicは2L+M、f=a²bである。R=3Lではない。
baseのfree divisor hxCがessential三であることと、S上のPoisson zero cubicの本質変数を混同しない。

全normal weight零symbolをcanonical base-Jacobian gaugeへ直すと、両源は同じ既存P0へ入る：

\[
 (P_{0,ac},P_{0,bc},P_{0,uc})=(-a²,2ab,2au),\qquad
 \text{他の全weight零括弧}=0.                     \tag{P0}
\]

nodeではE∧(2a∂c)を加える。
cuspではE∧(2a∂c)を引き、tensorを−1倍する。
いずれもprojective classの正当なEuler gauge／非零尺度であり、新たなparameterを反転しない。

## 2. node-flexの全higherを保持した表示

全normal-form磁気をR=Iu²+Juaとする。元六magneticへの正則inverseは前noteを使う。
(P0)を含む全tensorは

\[
\begin{array}{ll}
 P_{ab}=a(3v-u),&P_{ac}=-a²+2uv-2v²,\\
 P_{bc}=2ab+(u-v)c+Iu²+Jua,&P_{uc}=2au,\\
 P_{vb}=2v(v-u),&P_{vc}=P_{uv}=P_{au}=P_{av}=P_{bu}=0.
\end{array}
\]

raw元sourceの二次anchorと全magneticを単にrenameしてEulerを加えた式で、pointを固定していない。
P1にはJ、FULL9の中のP2にはIが入る。
二係数を一律にhigherのhomogeneous kernelへ押し込めず、inhomogeneous forcingを保持する。

## 3. cusp-flexの全higher

同じ(P0)に対し全非零括弧は

\[
 P_{ab}=-3av,\quad P_{ac}=-a²+3(u+v)²,
\]
\[
 P_{bc}=2ab+vc-Iu²-Jua,\quad P_{uc}=2au,
 \quad P_{vb}=-2v(v+u).
\]

他は零。nodeと同じweight零源だが、higherもdivisor配置も違う。
この二つを一actual orbitと同一視しない。元sourceの固定点列を増やしていない。

## 4. 旧constant-root分類との全仮定比較

`PLANE_NONREDUCED_M3_DOUBLE_SIMPLE_CONSTANT_ROOTS.md` §2の原全P0は

\[
 P_{ac}=-a²,\ P_{bc}=2ab,\ P_{cz}=(4-3e)az,
\]

及び全q connectionである。normal(z,w)=(u,v)と置けば(P0)は正確に
**e=2、alpha=gamma=delta=0**、すなわちLOG e2 whole lineの既存alpha零fiber。
G3はa²u²vというline全体を残すため、nodeとcuspの非零定数係数の差を捨てない。

旧全m三E定理は「actual D被約既約」を仮定する。
そこではfull Ore moduleのproper actualD因子を理由に多数の源を排除し、essential三分岐は**既約**freequintic八曲線へ送った。
本稿のDはhxCであり、そのproper factor枝を空／既八成分へ戻すことは出来ない。
同じ旧P0／full higher sourceに入るという事実だけではknown-parent positiveにならない。
逆に、必要な新数学はこの既存whole alpha零fiberとそのm二leading非零への変形に限定でき、
新P0点又は新linear Jacobi operatorを作る必要は無い。

## 5. cusp接線・割線のgeneric磁気はm四である

第三源は公開全normal form

\[
 P_{xz}=-hx,\ P_{xw}=3hx-3x²,\
 P_{yz}=3(h-x)²-3hy,\ P_{yw}=-9xy,
\]
\[
 P_{zw}=3xz+3(Jh²+Iy²)
\]
を使う。同じS、(a,b,c,u,v)で

\[
 P_{ab}=3(u-v)²-3ua,\quad P_{ac}=-9va,\quad
 P_{bc}=3vb+3Ju²+3Ia²,
\]
\[
 P_{vb}=-uv,\quad P_{vc}=3uv-3v²,
\]
他は零である。I≠0のS上zero cubicは3L、f=Ia³。
weight零はP0=3Ia²∂b∧∂c、normal columnsは全零。
actual Dは非零定数×uv[u²a−(u−v)³]なので、S沿い最小normal位数は**四**、先頭quarticはa u³v。
従ってm三constant00へ第三源を暗黙に含めない。I0はこのrank二plane入口の閉補集合として元A6へ残す。

旧m四31/F0,F1の全空性もactual D既約を必要とし、binary quinticのu/v因子を先に除いて正常形を作る。
本源はその除外されたbinary-factor境界なので、そのrank二空性を流用しない。
旧不可約G入口の成功を、現在のreduced1+1+3源のE完了と取り違えないための明示scopeである。

## 6. 次の最弱成功条件

最初の二源の既保存点はfullD101、射影接空間24。
同じ既点がnamed24成分T36/T37/T38のどれかへ入ることを実PGL曲線で示せれば、
局所次元下界24と接次元24で正則唯一成分になり、元全A6源の既約性から全I,Jをその成分へ戻せる。
新point列は不要。第三源の旧接次元29にはこの一anchor正則性橋を使えず、全lambda又はwhole sourceの曲線が要る。

MPS親への陽性が無い段階では、四次・五次normal jet又は別frameの点探索へ移らない。
旧sourceに戻した新二全表示と全磁気正規化の別作者手読解、
あるいはroot承認後のONE全150／Pf／実逆恒等式readbackだけを次の有限検証とする。
