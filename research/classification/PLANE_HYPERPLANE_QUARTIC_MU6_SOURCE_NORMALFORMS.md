# μ6/r2全十一源：全六磁気を一つの23軌道閉包へ戻す

2026-11-07、PRIMARY E。producer手橋とexact finite source batch、独立確認前。
全marked geometryは `PLANE_HYPERPLANE_QUARTIC_MU6_MARKED_LIST.md` による。
全source結果は `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_SOURCES_20261107.json`。
source一つごとのambient成分性・既知親membershipはまだ含めない。

## 全source batchで確認したこと

全十一候補と二E6 controlを一回の有限jobで確認した。
十一全てで、全28×18 gradient-syzygy行列がrank十六、全kernel二次元。
全二次Saito basisとactual markedH→h GL3及び逆を保存した。
その全基底からh-Casimir anchor二本を作り、全六磁気を保つ実sourceを構成した。
strict十Jacobi、projective三百五十Jacobi、実Pfaffianを全係数で確認した。
座標順(h,q1,y,q2,w)ではaffine Pfaffianはdet(U,V)であり、
(q1,q2,y,w)の順ではminusとなる。このorientation差はsourceやrankを変えない。

十一全てのcanonical cofactorは5×5 rank二・平方零で、hは双対像の外。
affine linear partはrank二、constantはその像内であり、正確にnil22-image/type6である。
二controlsも保持した。E6 a0 transverseは全quadratic anchor pairのdeterminant零、
非零a transverseは全六磁気でrank一nil211/translation零となる。

全十一で磁気shearの6×6作用行列Tはrank五。
この行列、小行列、全cokernel basisを公開reportに保存した。
計算は3.326秒、最大RSS73,456KiB。旧runnerの全FunctionDef ASTを保つisolated modeを使用した。
このPASSだけからcomponent数を増やさない。

## 全磁気のactual有限shear

任意のsourceを

```
P=partial_q1 wedge U+partial_q2 wedge V
 +(alpha*q1+beta*q2+R)*partial_q1 wedge partial_q2,
R=sum_j m_j*mon_j, mon=(h²,hy,hw,y²,yw,w²)
```

とする。q'=q+L、Lはbaseの全linear二成分、という真の可逆PGL shearの磁気作用は

```
delta R=-alpha*L1-beta*L2+U(L2)-V(L1)=T*l.
```

全Tはrank五、cokernelをI(m)=v*mとする。v_j!=0となる最初のjを選ぶ。
`m_normal=(I(m)/v_j)*e_j` とすればm_normal−mはTのimageに入る。
reportの非零5×5minor T[R,C]を使い

```
l_C=T[R,C]^-1*(m_normal-m)_R, l_other=0
```

と置けば、全六係数に対してT*l=m_normal−mである。
固定有理定数だけを反転した線形式なので、未知magnetic coefficientは一つも反転しない。
元の全sourceを正則・可逆なshearでR=(I/v_j)*mon_jへ送り、I=0も同時に含む。

I!=0ではq_new=q/(I/v_j)、tensor multiply I/v_jによりR=mon_jへ移る。
I=0は正則な実族R=t*mon_jのt=0端点。
従って各**全六mag source**のPGL飽和閉包は、選んだcanonical R=mon_jの一軌道閉包に等しい。
この逆変換は全magnetic quotientを保ち、単一係数切片の計算で全sourceを代用していない。

canonical Rは次の通り。

| 配置 | R |
|---|---|
| A6 cusp-secant bitangent | yw |
| A6 smooth flex | y² |
| A4+A2 cusps join | yw |
| A4+A2 A2 tangent | y² |
| A4+A2 A2-secant/smooth-tangent | y² |
| A4+A2 smooth flex | y² |
| 3A2 cusps join | yw |
| 3A2 cusp tangent | y² |
| 3A2 smooth bitangent | y² |
| E6 smooth flex | y² |
| E6 smooth bitangent | y² |

## 各canonical軌道の次元は正確に23

全underlying quarticのprojective automorphism groupは有限である。
A4+A2では異なる二cuspを固定し、normalizationのr-scaleを第三sectionのlinear/constant比が固定する。
A6ではnormalizationの同じ比較からr→zeta*r,zeta³=1だけ。
3A2では三cuspsの逆像を置換するS3だけ。
E6非零aではcuspを固定するP1変換のtranslationは第三sectionのst³項を生むため零、
remaining scaleはaを保つ条件で有限となる。
marked pairのautomorphismはその部分群なので有限である。

各実D=hFのdirectional-derivative kernelはintrinsicなfiber二平面。
全PGL5 stabilizerはこれと唯一のHを保つ。
infinitesimal base quotient actionは上の有限marked automorphismにより零で、
残りはq'=Aq+Lだけである。
projective tensor equality G_*P=sigma*PのanchorがA*M=sigma*Mを与え、A=sigma*Iとなる。
nonzero magnetic cokernelの比較はsigma²=sigmaなのでsigma=1。
残るtranslation kernelはdim kerT=1。
従ってcanonical pointのprojective stabilizerは正確に一次元、軌道は23次元。
各全source閉包も正確に23次元である。

## 正当化された次の唯一有限batch

各全sourceが一つのcanonical23軌道閉包であるため、上の十一点で一括して全150/350ambientを調べればよい。
D102ならsmooth23、D101なら唯一normalの非零二次余核が最弱で十分。
二次余核零又は複数normalならそのままfrontierとして保存し、深いjet・別anchorへ進めない。
D5+tangentの既済D101/Qcoker零点はこのbatchで再計算しない。
全既知111との相違及び独立review前に新componentを数えない。
