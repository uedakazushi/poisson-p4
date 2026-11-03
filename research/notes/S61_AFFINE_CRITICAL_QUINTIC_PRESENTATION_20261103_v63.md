# original critical Eのrank5有限表示と五次norm消去

2026年11月3日。一般affine fiberの全正次元支持を含むoriginal critical quotientを、既存のcyclic/quintic開で、**16方向座標上の5×680有限presentation**に置き換える。全一般fiber零次元性や五点尽くしは結論しない。

## 対象の開集合

元134式をFi(beta,q)、元45個のqについてのJacobianをJqとし、既存整数行列CによりE=det(C*Jq)とする。座標変更後のJacobianに置き換えない。actual quotientをA、critical quotientをAcrit=A/(E)とする。全正次元成分の被約支持はAcritに入るが、特異孤立点も含まれ得る。

以下で反転するのは、既存setup minors、方向chartのY1、cyclic determinant、monic quinticを作る2051 block determinant、free5代数のtのnormである。**gの判別式、normal-recognition行列Mのnorm、Jacobian normは不要**なので、旧H開の一部の補集合も残す。cyclic/high-block零集合を覆う主張ではない。

chartはY1=v1を取る。v1は凍結点で−3のp単元で、残る16方向はYrest−vrestである。旧Y1=1 chartの定数rescalingで、中心はY=v,t=1となる。Norm(t)を反転するとtによるdehomogenizationが正確にできる。Norm(t)=0でも一部のrootだけtが単元の場合があるが、それは本稿の開の外に保持する。

## 全critical quotientの小presentation

Bを五つのadmissible beta parameterと16方向を含む係数環とする。既存の有限inverse recipeで得るmonic quintic gに対しF=B[s]/gはfree rank5で、squarefree性は不要である。χ1031と元の全134 homogeneous Fiの余りをr0,...,r134とすれば、実affine代数はF/(r0,...,r134)である。χも残すため、free1031記号の余分な固有値を実根と同一視しない。

F内でtのinverseをadjugate/normにより作り、元のaffine45座標qaff(s)を復元する。その後に元q-Jacobianを評価して

    e(s)=rem_g E(beta,qaff(s))

を作る。これは既存有限inverseと45次determinantによる具体的な有限回路である。一般係数の巨大な展開表を保存したという主張ではない。

Acrit=F/(r0,...,r134,e)であるから、基底1,s,...,s4について全rj*s^kとe*s^k（0≤k≤4）を列とするMcritは

    B^680 → B^5 → Acrit → 0

という**全有限presentation**となる。free5代数の任意のideal倍はこれらの列でB-linearに生成され、全134式を保持する。supportはFitt0、すなわち全5×5 minorsのidealで正確に与えられる。そのunit検査は未実行である。

## Norm(e)による方向座標への消去

eの5×5乗算行列MeはMcritの最後の五列である。NE=det(Me)と置くと、Me*adj(Me)*1=NE*1はF内でNE=e*hを与える。従ってNEはAcritの正確なannihilatorで、critical支持の射影は方向変数だけのNE=0に入る。

NEは全free5のnormなので、actual quotientで消える余分なquintic rootも零集合を大きくし得る。従ってNE=0は必要条件でありsupport等号ではない。正確なsupportは全5×680 presentationで決まる。

## 新しいcenter certificate

凍結toric betaとY=vで、法32003のgは

    s^5+27940s^4+12935s^3+12457s²+3355s

で、相異なる根は0,1865,19521,6706,7974である。保存χのg可除性と、全29 coordinate polynomialの五根での値を再検査した。元45座標へ戻して全134式が0であることも直接確認した。

元q-Jacobianから計算したEの五値は18384,6387,6706,13504,28620で、その積は**NE(v)=9485≠0**である。CRTでdegree<5のe5を再構成し、Euclidean Bezout U*e5+V*g=1を全係数で確認した。Meとそのinverseの全25 entryも照合した。gがsplit squarefreeでdegree<5なので、これはfree5代数の全係数identityである。

中心では元675 residual列が全て0、最後のMe五列が可逆なので、critical quotientは厳密に0である。これは一般NEの非零性と中心を含むnorm開を与えるが、全方向でunitとは言わない。

## 証明書と残余

* `code/s61_affine_critical_quintic_presentation_20261103_v63.py`
* `reports/S61_AFFINE_CRITICAL_QUINTIC_PRESENTATION_20261103_v63/{report.json,finite-recipe.json,center-certificate.json}`

新規検査1.09秒。一般presentationはexact finite-inverse recipe、数値certificateは凍結F32003での全係数証人である。quintic開内の未知critical支持は5×680 Fitting問題として残る。cyclic/high-block/t-norm補集合、全一般affine零次元性、未知孤立点、全五点尽くしは未解決である。
