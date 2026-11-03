# 基底と個々の被約成分の pinching 完成環：独立監査

2026年11月3日。基底 Z と CA/CB/N の二つの frozen producer を新しい作業ディレクトリで実行し、
全 QQ 証明書が実行秒数以外一致した。個々の全61座標・全134式・150係数・projective Jacobi
の producer 再実行を含む。さらに独立コードで有限二形式作用、全 toric 消去、実 β 座標、
全45個の固定 target K における一次奇数係数と ideal を検査した。

受理した環は基底で C[[u7,s,p²,p³,sp]]、個々の被約成分で

    Ac=C[[u10,s1,…,sc,p²,p³,s1p,…,scp]],
    N:A7, CA:A8, CB:A8.

全完成環への論理を特に確認した。proper birational な滑らかな旗入射は基底旗の有限長2
ファイバーの近傍で有限正規化となる。その完成への単射と、同変な H7×slice の厳密な積分解を
使うので、任意の線形 section と scheme image を交換する議論ではない。
source の二変数族が全横断正規化を形式的に parametrise することは、七つの固定 β 座標の微分が
階数7で、その核を ds,dp が張ることから従う。

実際の target 座標 K,s を regular 座標とすると B=R[[p]]、R=C[[s,K−K(Q)]] である。
A は R,p²,p³ を含む。D=R[[p²]]⊂A⊂B とおくと B は D 上有限自由階数2で、A は D 部分加群
として有限かつ閉である。従って p²B⊂A は全形式的無限和について成立し、
全像は A=R+Ip+p²B と厳密に決まる。一次奇数係数の検査を全環の結論へ上げる根拠はこれであり、
有限 jet の一致だけを使っていない。三族の選択 q block の determinant は p,s に依存しない
非零定数であることも独立に確認した。

I は N で (s,k1,k3,k4,k6,k7,k8)、CA/CB ではこれに r を加えた ideal である。
全61座標を p=s=r=k1=k3=k4=k6=k7=k8=0 に制限すると三族は文字通り一致する。
共通被約非正規層は自由 k0,k2,k5 と群方向7の滑らかな10次元層である。
正規化商 C[[u10]]p、conductor (I,p²)、depth11、multiplicity2、embedding dimension26/28
は完全列と有限 birational multiplicity の加法性から従う。

指定点は a=1 の Q であり、a=2/9 の共鳴点ではない。基底の torus 正規化は C 上の非零五乗根を
使う正確な同変移送である。指定 Q の長さ2正規化ファイバーから I が全て零になることが分かり、
上の非正規 model が適用される。

範囲は個々の被約 S61 成分の局所環である。全 S61 の局所環、交叉の全 scheme ideal、
他の枝をこの結果から決めない。元の射影模型への移送は generic Poisson rank4 開部分であり、
rank≤2 境界の全域 R1 や交叉成分の不存在は未検査である。
明示的な範囲追記 AMBIENT_INTERSECTION_SYMPLECTIC_SCOPE_20261103_v63.md も監査入力に含めた。

再現コードは code/ambient_pinched_models_review_20261103_v63.py、結果は
reports/AMBIENT_PINCHED_MODELS_20261103_v63_REVIEW.json に保存する。
