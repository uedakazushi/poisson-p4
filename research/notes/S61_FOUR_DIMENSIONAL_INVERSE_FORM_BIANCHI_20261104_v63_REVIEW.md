# 4D inverse-form/Bianchi の独立検証

2026-11-04。主producerはimportも再実行もしない。独自のsparse integer/Fraction polynomial辞書を使い、110個の独立0/1/対称2jetとformal h,dhについて、全係数で検証した。主report/source/certificateのhashと、保存されたP,A,p,dp,X,j,E-dual,三operatorの全式も独自AST parserで照合する。

## 符号と全普遍恒等式

PfaffianはPを2formとしてP∧P/2から、Aはnegative Hodge complementから構成した。外微分はtarget multi-indexごとの交代和、wedgeはtargetのsubset shuffleとして作った。主コードのloop順・adjugateのhardcoded matrix・SymPy polynomial ringは使わない。

Schouten符号は独立なGrassmann計算で固定した。左θ微分を用いるbivector Pについて

    [P,H]=−sum_i (partial_theta_i P wedge partial_x_i H
                              +partial_x_i P wedge partial_theta_i H)

を計算する。これにより[P,P]=2Jと、Hのdualをhにした[P,H]の4vector係数がX·h−P_ij partial_i h_jになることを双方確認した。

PA=AP=pI、j=−dp−AX、dual(E)=−Pj、dE=2dp∧dA、p dE−2dp∧E=0、全secondjetのactual Bianchi、formal全110変数の

    D_P(h)=p B_P(h)−2j^T P h

を全係数で確認した。特にdEの符号はpositive2であり、dual3formから4formへの外微分はdivergenceとなるので、主記述のcleared operatorと一致する。

## 元source generatorへのtransport

j=QFと書き、S=Q^T P QとするとSは交代行列である。差cofactorのb成分は

    (−2F^T S)_b = sum_(a<b) −2S_ab F_a
                               +sum_(b<c)2S_bc F_c。

従ってcofactor vector全体は

    sum_(a<b) −2S_ab (F_a e_b−F_b e_a)

であり、任意個数Nの標準Koszulの多項式結合となる。Qの正方性・可逆性は使わない。Qが空間変数に依存してもoperator差のLeibniz項は既に同じ普遍式に含まれ、Koszul transport自体はpointwise polynomial identityである。独自コードは別のformal4×3 QとF3でも全cofactorを確認するが、任意Nの証明は上の添字式である。

空間係数を取る操作は単なる有限convolutionなので、S[x]上のBianchi/Koszul combinationは各空間係数でも対応するsource moduleに入る。full134/135 correspondenceは既存の独立受理済みfull operator監査を依存として用い、今回巨大展開を再実行しない。

## 非単項式saturationのscope

block diagonal P01=a,P23=fにすると、差のh0係数は実際に−2a² partial_1 f。独立jet変数fのfactorを持たないので、p=afでは割れない。したがってKoszulを除く前のcofactor vectorをpで割る操作は一般にはpolynomial source vectorを与えない。Koszulを除けば残りはp Bであり、この特定操作からprimitiveに戻るものはBそのものである。

これは一般のmodule saturationを計算した結果ではない。localizationへのpreimageはM:g^infinityであり、SがdomainならM[g^-1]∩S^Nと記述できる。新candidateには全cofactorのpolynomial性とg^k v∈Mが必要で、既存moduleに新しくないと断定するには別途v∈M、真に新しいと言うにはv∉Mが必要である。空間p(x)を逆にするなら最初にS[x]で行う。p(x)をsource係数環の単元として扱わない。

従って主packetの限定的method gateを受理する。一般nonmonomial saturation、全kernel、残る43方向、full80、四次239、一般beta尽くしのnegative/positive結論には使わない。主source bytesは変更しない。
