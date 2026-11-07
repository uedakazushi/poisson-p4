# rank零plane・globalV二parameter四族：三族全閉鎖とnil線残余

2026-11-07、主義務G。producer exact whole-source証明。
全normalizationとinverseはnormal半単純/冪零first-jet notesによる。
以下はいずれもfixed-plane affine symbol base十一で、全Π1三十・Π2九を保持する。

| whole family | T | V | K |
|---|---|---|---|
| normalss MIXH(u,k) | uH+N | H | kzw |
| normalss MIXN(u,k) | H+uN | N | kzw |
| normalss NILREG(u,k) | uJ+J² | J | kzw |
| normalnil b零MIX(r,k) | H+rN | H+N | kz²、r≠1 |

H=diag(1,1,−2)、N=E01、J=E01+E12は係数行列であり、Eijはfield xj∂xi。
u,r,kを固定標本へ置き換えず、full QQ[tau,eta] sourceで全三百五十coefficientsを使った。
generic higher kernel upper dにはambient16+d、whole wallでは15+d、fixed pointでは14+dを適用する。
毎whole lineをQQ(t)のkernel/minorと全exact finite rootsで覆い、line intersectionも各lineのgeneric patch又はrootに入る。

## 1. 全MIXH/MIXNと全normalnil mixed

`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TWOPARAMETER_20261107.json`に全operators、
complete kernels/minors、全maximal-minor factorsとwholeline/rootを保持する。
normalss MIXHは全十linear wallsと全pointを含めimage dimension二十三以下。
normalss MIXNは全九linear wallsと全pointを含め二十二以下。
generic kernelsは両者とも六/零である。

normalnil b零MIXの三linear wallsも全pointを含め二十二以下。
唯一nonlinear polynomial factorは

`eta*tau−tau+3=0`。

このwhole curve自身がtau≠0を強制し、全体を`eta=1−3/tau`でparametrize出来る。
全D1/D2/shear matrixをこの式へ引き戻すと、分母はtauだけ。
三つのfull matrixをnonzero scalar tauで掛けることで、rank/kernelを変えずQQ[tau]上へ戻した。
これはsourceの不明な補集合を削るlocalizationではなく、tau零のcurve部分は空である。
全minor rootsはtau=0,1だけ。前者は空、後者はT,V dependentでm三以上の別義務。
それ以外のwhole curveではkernel六/零、source base十より全像二十一以下。

exact証明書は`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_HYPERBOLA_20261107.json`。
したがって全normalnil b零mixed m二sourceも、低parameter noteのspecial T/V行を併せてG非優越である。
全normalnil b非零は既に別packetで二十二以下。
全normalnil b零regular-nil/common-kernel m二sourceも低parameter closureで二十四以下。
normalnil b零Cartan三parameterだけは本三つの結論へ含めない。

## 2. normalss NILREGの正確な残余

k非零では全四linear walls u=±1,±3のnonzero-k点を含め、上のlinear-kernel methodで閉じる。
k零のwhole u-lineではgeneric d=11、image upper二十六を残す。
全minor exceptional valuesはu=±1/3,±1,±3である。
u=±1/3のfixed-point linear upperは二十五、u=±1,±3では二十六。
これらは真のPoisson-family dimensionを証明した値ではない。

全actual shearは六で、generic weight一kernel八から取ったtransverse sliceは二変数。
weight二kernel三を全て戻すとgeneric normalized sourceは五変数であり、
compatibility binary quadricsのheight二でG≤24に十分となる。
固定u=±1,±3ではweight一slice三変数にheight二、u=±1/3では二変数にheight一が十分。
全QQ(u) sourceのbasis/minor例外も復元する有限binary-resultant/point proofが次の義務。
未実行のheightを仮定してこのNILREG全源を閉じたとは主張しない。

## 3. custodyと残るwhole sources

二parameter実行は240秒/1536MiBのfinite contract、実測19.005秒/79,004KiB。
hyperbola全source実行は90秒/1536MiB、実測2.519秒/65,040KiB。
reusable runnerへisolated modesをappendし、各回全旧FunctionDef ASTを直前immutable reportのsourceと比較した。
旧二十三file manifestの可変九入力（runner84f1＋八jobs）は
`inputs/plane_rank_zero_global_lowparameter_frozen_inputs_20261107.json`へhash付きで保存済み。
既存certificateのbytesをcurrent runnerへ書き直していない。

normalss generic MIX三parameter、normalnil b零Cartan三parameter、上記normalss NILREGのK零線、
rank零scalar co-Higgs、pointwise-dependent T,V、高mのO(1)rank一及び全leading-Pf零復活sourceは別残余。
本noteは全k零/globalV、全m三以上又は全Gの尽くしではない。
