# K4a の actual24 親は T36/T38 でない：全order二零点の空性

2026-11-08、PRIMARY E。全大域 operator の producer PASS、別作者 operator 再生は未完。
入力と intrinsic 条件の全証明は `K4A_ORDER_TWO_ZERO_LOCUS.md`。
全非零eが一つのactual24orbitに入る、独立受理済みliteral curveの整数parameter e=1 を使った。
新しいrank点、q chart、normal又はprofileを選んでいない。

## 全 P4 上の scheme を一括排除

全十 cubic `T=E∧P` の全五 first derivatives、計50 quadraticを残した。
15個のdegree二monomialに対する50×15の係数行列は、QQ上rank15、kernel零であった（0.335秒）。
全50式・全行列・15独立行列・非零QQ determinantを
`CLASSIFICATION_K4A_ORDER_TWO_ZERO_LOCUS_20261108.json` に保存した。

従ってC上でもこの行列はinjective。
もしq∈P4がorder≥2零点なら、その非零monomial vector `(q_i*q_j)` がkernelに入る。
少なくとも一つのq_i²は非零なのでこれは不可能。
すなわちactual dense24orbitにorder≥2零点は一つもない。
有限点標本、実点又はQQ点だけの空性ではない。

## Completed T36/T38 の全閉包に必要な性質

両受理completed sourceの全括弧はx0に依存せず、斉次二次係数である。
q=[1:0:0:0:0]で全括弧とその一次微分が零だから、
`∂(E∧P)(q)=0`、従ってprojectiveΠのintrinsic order二零点である。
一般period/smoothing parameterを一つも固定する必要はない。
period分母が非零の全generic sourceに成立し、全PGL飽和にも成立する。

全tensorとqのincidence `∂(E∧P)(q)=0` は閉条件で、qのparameterはproperなP4。
従って「order二零点を持つ」は全tensor空間の閉条件である。
accepted fullfamily component theoremにより、T36/T38の全component closureでもこの性質が成立する。
ここでgenericのknownpointから別sourceへ帰属を外挿していない。

actual K4a親のdense24orbitはこの必要閉条件を破るから、親成分はT36/T38と異なる。
canonical K4a P0は独立確認済みsmooth unique24点なので、もしP0がT36/T38に属せば
unique componentはその成分となり、denseorbitも同じ必要条件を満たして矛盾する。
従って元の全一orbit K4a source自身もT36/T38の外にある。
sourceの境界全点までの非所属は主張しない。

T37又は未登録24親は本条件で除外しない。
他の既知24親を閉因子/本質変数条件で除く監査は、統合担当の全prior比較に依存する別義務。
新成分の計数は零のまま維持し、次はT37のactual source認識又は別の全相違証明へ進む。
