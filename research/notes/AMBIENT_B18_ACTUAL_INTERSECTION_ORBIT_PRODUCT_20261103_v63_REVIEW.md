# B18 actual Dのorbit productとdepth範囲：独立概念監査

2026-11-03。対象は指定P0のactual intersection scheme D=CA∩CB。
既受理のdim28/edim39/e2/depth≤20を入力とし、十九smooth factorによる下限だけを加える。
exact depthや残るsliceのsocleは決定しない。

## nilpotentを保持するscheme product

CAとCBはそれぞれGL不変な実被約成分なので、その**actual ideal**はGL不変。
和idealで定義されるDもschemeとしてGL不変であり、D_redだけの不変性を使っていない。
canonical126係数空間はsmoothかつGL安定である。
P0のorbit tangentの十九独立列に対応するformal group section H19と、十九独立linear coordinateを
固定する補affine sliceを選ぶ。sectionはsubgroupである必要がない。
`(g,s)->H(g)*s` のambient tangent mapはisomorphism。
formal inverse function theoremにより、ambient completionはgroup19×slice107へ移る。

G不変idealをこの作用mapへ引き戻すと、元idealはgroup parameterを含まないslice idealと正確に一致する。
逆向きにもH(g)^-1の作用を使えるため、単なる点集合の一致でなくidealの両包含が成立する。
従って全nilpotent・embedded moduleを保ったまま

    D_hat ≅ D_slice[[g1,...,g19]]

となる。source imageと任意base changeの交換や、Dのreduced性・CM性を仮定しない。

## 既存物理行列からの独立有限検査

元P0の四つの非零二次係数から、外積indexとinverse-coordinate monomialの微分で全25 GL列を再構成する。
producerのLie actionは再実行しない。十九minorの非零と全25列のそのspanへの包含を完全行列等式で確認する。
`E_yy` の列がP0自身に一致することも確認するので、radial directionはこのorbit factorに含まれる。
元の150×39 actual tangent監査を依存として固定し、その全入力hashを照合する。

## depth・射影化と未決部分

sliceはdim9、edim20、multiplicity2。
十九power-series variableはregular sequenceなのでdepth D=19+depth D_slice≥19。
incomparableなminimal R20 componentからの既受理上限depth D≤20と合わせ、

    depth D ∈ {19,20},   depth D_slice ∈ {0,1}。

embedded vertexあるいはslice socleを排除していないため、どちらか一方へ決めていない。
新しいfull local ring、D_red tangent、nilpotentのlinear visibilityも推論しない。

P0は非零cone点でq106=1を持つ。
Gmのscalar作用はそこで自由で、q106によるdehomogenizationはregular factor一つを除く。
必要ならorbit basisをE_yyと他十八列へ取り直せるため、projectiveではsmooth factor18と同じsliceとなる。
従ってprojectiveのdim27/edim38/e2はそのまま、depth範囲は18〜19。
これはZariski近傍の全点・他のB18境界へexact depthを拡張する主張ではない。

独立結果はflat report
`reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json` に保存する。
