# 全coreの高rank比例case排除：独立監査

2026-11-03。author packet `S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63`を
同じ実四quadricから独立に再構成して受理した。
結果はa*P*J≠0での比例排除であり、全raw rank24/product41の証明ではない。

## 実sourceとglobal rankの確認

凍結archiveと全input hashを照合し、804個の実quadric係数をQQ上で再構成した。
W=x5=0でのQ1=A*B/6075、Q2-combination=40*x0*B/6075と
3Sc−7Vc=20Pを文字通り確認した。

split座標(r,u,v,t,R)から元の五座標へ戻すlinear matrixのdetは**1/(2a)**。
従ってa≠0で実可逆であり、4Q0=uR−5vtは単なる非可逆pullbackの恒等式ではない。
W=0上でQ0はrank4で、radicalはr方向。
P≠0ならSc,Vcが同時に零にはならず、Bは非零なhyperplaneである。
一つのhyperplaneへの制限でrankは高々2しか落ちないため、Q0|B=W=0のrankは全coreで少なくとも2。
rank0/1の場合をこの議論の外へ落としていない。

## DB≠0の実係数unit

DB=20P²+3Pe+20P+180a、br=−20a*DBとする。
a²Q3の実5×5二次係数行列Hは、authorのHessianを呼ばず、二次式の各係数から直接再構成した。
全25entryとBの全五係数は保存されたQQ多項式に一致した。

独立に確認した恒等式は

    Hrr*bU²−2br*HrU*bU+br²*HUU = −a^5*br²*J/12.

DB≠0ならrをB=0で除去でき、元の**unscaled Q3**のu²係数は−a³J/12。
a,Jの既知単元だけで非零となる。
同じ制限Q0にはu²項がないため、Q3はQ0のscalar倍ではない。scalar=0も含めて排除する。
このbranchでL=20P+3e+20をinvertしていない。

## DB=0の全cross-minorとQQ Bezout

DB=0ではa=−P*L/180。
a,Pが単元なので、ここで初めてL=−180a/Pが既知単元となる。
Hrr=0を全QQ[P,e]で確認した。
Q0はrを含まないので、Q3|BがQ0|Bに比例するにはr-linear row hがBのmain row bに比例し、
全cross-minor h_i*b_j−h_j*b_iが零となることが必要である。

保存された五つのGiについて、実際のcross-minorから

    cross_minor_i = P^k_i * L^l_i * J^m_i * Gi / scalar_i

となる全係数恒等式を確認した。
記録されたpowerは全て非負で、除く因子はP,L,Jと非零QQ定数だけ。
未知factor、DB、N、minor自身などを勝手にinvertしていない。

さらに、保存されたcofactorを独自のFraction sparse多項式積で掛け合わせ、

    sum Ci(P,e)*Gi(P,e)=1

をQQ[P,e]内で文字通り確認した。
従って全五minorが同時に零になる点は存在せず、DB=0でも比例は不可能である。
標本、有限体、GB、CRTは用いていない。

## gcd-gateへの正確な帰結

複素数上でrank3/4の二次式はirreducible。
同次数の別の二次式と非定数factorを共有すれば、両者は比例する。
それは上で全a*P*J openから排除された。

残るrank2 common-factor例外は、既存の実dual-norm/必要minor計算により有限閉cover Sigmaへ入る。
別の独立監査済み有限V-model Nakayama証明は、Sigmaの全座標代数上でraw24×56のQQ全射性を示している。
従ってFULL U上の**全gcd-gate failureはraw-good**である。
これはgcd failureが存在しないという主張ではなく、その例外もraw rank障害にはならないという結論。

pair gateとfive-quadric rank gateは残る。
元product41、quartic、ambient removalや全core raw24を、このreviewから完了と扱わない。
有限Σに限定したmodp→QQのfinite-proper論証と、今回のglobal QQ比例排除を混同しない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/s61_nonsquare_gcd_proportional_review_20261103_v63.py --output <new-output-directory>

全再構成Hとcross-minor:
`reports/S61_NONSQUARE_GCD_GATE_PROPORTIONAL_EXCLUSION_20261103_v63_REVIEW/`。
実行約2秒、制限120秒・2.2GB。凍結author sourceやwork出力は変更していない。
