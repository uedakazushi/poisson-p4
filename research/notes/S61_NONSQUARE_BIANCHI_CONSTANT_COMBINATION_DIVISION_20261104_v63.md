# QQ定数結合と共通一変数因子によるweight4候補

2026-11-04。凍結済みhigher Bianchiのweight5/6/7/8群（20/22/18/11行）を使う。元operatorや135式を再構成せず、weight0 nuisance19/23だけをliteral Tnormal/U zeroの確認後に代入する。

群のtotalweightをw、factor候補local variableをv（weight(v)=w−4）とする。全135cofactorのうちvで割れない全monomialの係数を、a,e,dの指数ごとに**QQ係数へ分離**し、行結合の未知QQ定数≤22個に対するexact homogeneous linear systemを作る。数値parameter代入やQQ(a,e,d)への拡大は使わない。

同じ群の全cofactor mapのkernelも求める。これは恒等的に零になる行結合であり、divisibility kernelからこの部分をquotientする。保存するfree-coordinate kernel vectorsの全constraint積が零であることを検査し、選択constraint rowsが正確なrankを与える。divided imageに対する別のexact eliminationでnonzero quotient dimensionも照合する。

非零候補が得られた場合は全combined coefficientsがvで割れること、quotient totalweightが4であることを確認する。F5/F6を含む全135列をdivisibility constraintsに保持する。除算後のcofactorsと、同じnuisance代入をした全135実式とのliteral polynomial恒等式を確認してから保存する。tG6のtだけを個別にunitとして除かない。

候補の保存はsource-levelの結果であり、normal graph・T34/Tn frame transport・354→319への輸送・新pivot・独立方向数は未検査である。候補がない場合も、negativeは**同じweight群のQQ定数結合を一つのactual variableで割るclassだけ**。parameter polynomial係数の結合、general gcd、module saturation、不存在一般論を含まない。

実行上限60秒・2GiB、artifact32MiB。各weight群の完了checkpointを保存し、cap到達時には失敗状態を保存して自動再試行しない。GB、新しいfield/point、元operator再計算は行わない。
