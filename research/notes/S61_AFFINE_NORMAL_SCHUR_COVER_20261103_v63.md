# Generic six columnsによる未知残余の15個のscalar chart

2026-11-03。全未知アフィン残余Bを失わない有限coverを構成し、各chartで5-target moduleを一つのscalar idealへ移す。これは構造的な定理とfinite recipeであり、15個のunitそのものは計算していない。

## 1. 行列のrankと未知base

K=C(beta16)とし、既証明の非空beta開の共通部分で考える。formal R'=K[n42,S1,S2,H]上のglobal monic quinticは、実際の全AをPhi:R'^670→R'^5のcokernelとして表す。known-five splitとnormal閉埋込みより、R'-moduleとして

    A ≅ (R'/m)^5 ⊕ R'/I_B

であり、Bのsupportはmを避ける。ここでE5は一般にK^5という代数ではなく、K-dimension5のmoduleとして上式を使う。

従ってB上ではPhiのfiber rankはちょうど4であり、imageはrank4のlocally free moduleとなる。既存の実finite two-coordinate projectionよりdim B≤2である。これはgeneric betaの主張であり、凍結有限体の全special fiberに適用したとはしない。

quintic target basisは(1,h,h²,h³,h⁴)であり、B上でのquotient rowは(1,h,h²,h³,h⁴)である。R'→Bが全射なのでcoker(Phi tensor B)=Bで、このrowのkernelはbottom四座標への射影と同型なfree rank4 B-moduleである。したがって**Phiの固定bottom四行(1,2,3,4)はBの全点でrank4**である。五つのrow選択は不要となる。このB上の自由性をR'上のAの自由性と混同しない。

## 2. 六個の定数column combinations

670×6行列Lambdaのentriesを4020個の独立なconstantとし、faithful field extension L=K(lambda)を取る。Psi=Phi*Lambdaとする。

固定したBの幾何点ではPhiのimageは4次元なので、Psiがrank≤3となるLambdaの条件の余次元は

    (4−3)(6−3)=3

である。Bの次元は≤2だから、全悪いincidenceの次元はparameter空間より少なくとも一つ小さい。その像の閉包は真の閉集合であり、generic Lambdaは**Bの全点**でrank4を保つ。generic pointだけのrank検査ではない。非被約部分も、fiber生成とNakayamaにより失わない。

よってPsiの4×4minorsはB_Lを覆う。fixed bottom four rows、column選択15通りの計15個である。全15selectorを保存する。特定の数値QQ Lambdaが検証済みとは主張せず、明示的なtranscendental constant extensionを使う。このfield extensionはfaithfully flatなので、後のunit判定を元Kへ戻せる。

## 3. 各chartのexact scalar ideal

固定4row R=(1,2,3,4)と4column Cを選びdelta=det Psi_(R,C)を反転する。原点ではPhiの全entriesが零だからdeltaも零であり、このchartはknown-five factorを含まない。従ってA_delta=B_deltaである。

四本のpivot columnは元670列の定数結合なので、それらを追加してもmodule imageは変わらない。四つのtarget generatorを解く通常のSchur eliminationにより、残る一つのgeneratorのrelationsは

    g_j=det[Psi_(all5rows,C),Phi_(all5rows,j)], j=0..669

で与えられる。deltaは単元なので、これらのdeterminantsとSchur coefficientsは同じidealを生成する。従って

    B_L localizeddelta ≅ Spec R'_L[delta^-1]/(g_0,...,g_669)

がスキームとして成立する。生成数は≤670、parameter次数はdelta≤8、g_j≤10である。全134式のmodule columnsを保持するので、未知点やnilpotentsを省略していない。

これは一つのcyclic-openだけを取ってdeterminant零の補集合を捨てる方法ではない。15個が全未知supportを覆うことは第2節で証明済みである。なおPhiのrank4 locusはunknown、rank5 locusはsupport外であり、scalar equationsがその区別を保持する。

## 4. 正確な未完了課題

したがってBが空であることは、上の**全15個のlocalized scalar idealsがL上でunit**であることと同値である。いずれのunitも今回計算していない。有限field witnessによるgeneric unit推論も行わない。

このrecipeは新global5×670 moduleのexact support reductionであり、full42/45-variable GBを再試行した結果ではない。系の次数を上げる探索、minor集合の巨大展開、未知beta-finitenessの仮定は使わない。

保存: `code/s61_affine_normal_schur_cover_20261103_v63.py` と `reports/S61_AFFINE_NORMAL_SCHUR_COVER_20261103_v63/`。receiptは依存hash・数値次元・全selector・finite equation recipeを検査する。幾何学的証明は本noteである。
