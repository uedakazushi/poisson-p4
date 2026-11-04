# κ14と三raw交叉minorの一回検査

2026-11-04。凍結source gate b9de7b3aの合同73演算・75284係数枠・6.61MBに基づく。
四つのκ係数とanchor列は既存exact公開値を用い、current15四行[2,4,5,8]の列14のみを
元sourceから計算する。他列、追加行列式、点や素数のfamily検査は行わない。

κ14と三つのraw交叉minorを同じsourceから別々に計算し、全列14値・κ・anchor値と
四つの答えを保存する。source cache解放後、逆順四項積でκを再構成し、三raw交叉も
異なる積順で答えとの差が零になることを独立に検査する。QQ分母、a-Laurent offset、
全source支持上界と実出力byteを実行時に確認する。

κ14零だけならraw比例とは言わず、heavy二列を含むspan補正を全十五行へ適用して
初めて先の裸七borderを使える。三raw交叉が全て零なら、既存非零anchorと整域性により
raw列14もanchor列に比例し、既知六列と合わせて選択七pure列rank一となる。
他列への推論、残余三行rankや全FULLU cubic57は行わない。
