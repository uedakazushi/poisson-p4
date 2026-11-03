# Q_bのCB meeting double・actual D二次dataの独立監査

2026-11-03。author packet `AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63`を独立に再構成した。
全QQ(b)物理二次jet、rankと二次cokernel unitを受理する。
CBの二枝meetingの全型は受理するが、actual CA/CB intersection Dを先行P_aの型へラベル交換することは受理しない。

## 計算の独立性

全22個のCB enlarged canonical source basisを、3binary f、3binary g、16Hamiltonianのraw tensorから再構成し、保存済み全3300係数と照合した。
実fixed19、ratio、全flag9の微分は、この係数basisと独立の外積・二次monomial Lie作用から作った。
各29 source列を元の同一Q_bへ戻す際は、5×5定数行列の外積と逆二次座標の係数変換を使った。

本監査の二次jet評価器は、tの次数0,1,2を持つQQ(b)係数配列の有限convolutionである。
source部分は

    q0=canonical(rawCB(F,0,H)),
    q1=canonical(rawCB(deltaF,deltaMu*F,deltaH)),
    q2=canonical(rawCB(0,deltaMu*deltaF,0))

を直接含む。
exp(tA)とその逆の三係数を作り、行列積が恒等行列になることを確認した上で、
外積と逆二次monomial代入を係数配列だけで計算した。
authorの巨大な座標式展開やLie Hessian公式を呼び出していない。
両AはA^5=0で、実有限GL exponential curveを定義する。
今回の二方向ではdeltaMu=0だが、mixed項を零速度だからという理由で評価器から省略せず、非零mixed方向のsource係数も別に確認した。

保存した全150係数について0次・共通1次・真2次差が独立計算と一致する。
source方向の分母は非零有理定数とbの冪だけ。
全input hashを照合し、依存する独立監査receiptのPASSも確認した。
さらにauthor producerを新規directoryへ再実行し、elapsed time以外のcertificateとreportが完全一致することを確認した。

## 受理した厳密なrankとunit

CB両sourceとCA sourceはそれぞれrank29であり、QQ(b)上で次の値を得た。

    rank(CB1,CA)=36, rank(CB2,CA)=30,
    rank(CB1,CB2)=37,
    rank(CB1,CB2,L)=38,
    rank(CB1,CB2,L,CA)=38,
    rank(CA,v)=29, rank(CB1,CB2,CA)=37.

疎なfunctional e63（零始まり）は両CB接空間の全58列を恒等的に殺し、

    e63(L)=24*b^4

である。有限標本の非零性ではなく、全QQ(b)恒等式である。
方向の分母とこのunitにより、二次非零性は全b≠0へspreadする。

## なぜ全CB meetingがdouble S20になるか

既存の全旗監査により、全b≠0でCBは二つのsmooth29 branchを持つ。
実GL×b写像のrank20とcross-tangent rank37も、既存の単項式minorとGL-equivarianceにより全source近傍で一定である。
constant-rank theoremは両branchに入る実smooth20 germ Sを与える。

meetingの接次元はS上の全点で21なので、一次式を消去したsmooth21 ambientでSをrho=0とすると、meeting ideal Jは(rho²)に入る。
一点だけのrankから二重因子を推論したものではない。
独立に確認した24*b^4はnormal rho方向の二次式が単元係数を持つことを示し、J=(rho²)となる。
従って全b≠0で全meeting ringは

    C[[s1,...,s20,rho]]/(rho²).

全二branchの実imageはこのmeeting上のfiber productである。
先行の独立なkernel-completeness normal-form証明は環の抽象型だけを使用するため、そのまま適用できる。
共有s20を省く全81関係は

    y_i*z_j, (tau−rho²)*y_i, tau*z_j, tau*(tau−rho²)

であり、rho→(rho,rho), y_i→(y_i,0), z_j→(0,z_j), tau→(rho²,0)。
全81恒等式も本監査で確認した。
形式的normal form f(rho)+tau*g(rho)+Y(rho,y)+Z(rho,z)の一意性により全kernelとなる。
従ってcone dim29/edim38/depth21/e2、projective dim28/edim37/depth20/e2を受理する。
37はbranch接空間の和であり、whole CB edimではない。

## actual Dの非対称性と限定

CB2とCAはsmooth N28を含み、接交叉も28なので、その実intersectionはNである。
一方CB1とCAの接交叉は22であり、共通extra方向vもCA、従ってNの接空間内にある。
先行P_aで用いた、共通extra方向がNへ横断的という条件は成立しない。
よって旧full-D ring、nilmodule、長さ、primary decompositionのラベル交換は不可。

全CA接空間は既に二CB branch接空間の和へ入るため、actual Dのgeneric embedding dimensionは29。
既存の全交叉の被約台N∪R、R∩N=B18、B18⊂Sing(N)、およびこの族でのN smoothnessから、ここでの被約台はN28である。
従ってactual Dは一般にdim28/edim29で非被約となる。
このscopeを超えるactual Dの全局所型や二normalの高次obstructionは監査していない。
後者は別agentの継続課題であり、本監査で重複探索していない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CB_self_pair_actual_D_twojet_review_20261103_v63.py --output <new-output-directory>

独立物理jetと全rank行列:
`reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/`。
約9秒、制限120秒・2.2GB、GBなし。
