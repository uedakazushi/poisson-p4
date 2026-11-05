# HEAD29/30 の実cache合成：独立値検証

独立検証は PASS した。報告 SHA256 は
`bffdb3cb50ea194531b9e237da68e2612872aa0fbdc4edb32f7f9b2db8a5b058`。
primary `29133a93…`、SOURCE/code peer `985e6328…`、元compiled graphと全入力hashを固定した。
新しいimplicit解、physical/gauge jet、別headの値は求めていない。

独立実装は元producer・workerをimportしていない。
九変数指数をradix5整数へ符号化した疎な多項式dictionaryを使い、
積の前にtotal degreeを切ることでradix carryを避けた。
逆元は符号付きnilpotent幾何級数を作り、そのpowerの加算順も逆にした。
積・和の係数loopとpostの項順もproducerと逆順にした。
元n/d定数と固定105 lower入力から全8081 nodeを再計算した結果、
chosen lower-extensionの全node値hashが一致した。
全逆元の左右積も確認した。

delta用の二rootは内部63精度で計算する。
共有bB63 root38648もその祖先なので63精度を持つが、
bBとしての出力はcapture時にplane10へ射影した。
全safe delta63係数と八つのtrusted bB plane10配列を、独立計算終了後にprimaryと比較し、
全て一致した。既存C2、pureell、delta四次二係数、bB63のr ell²係数とのoverlapも確認した。

さらに、独立delta出力の全14個のzi*r/zi*ell二次係数を、
受理済みC2のliteral nonzz係数と比較し、全て厳密に零だった。
特にdelta63_z1ell=0を直接確認したので、HEAD30 generic convolutionで生じた
`−827/63360*delta63_z1ell` の一項も実値で消える。
これは形式変数を根拠なく零に置く処理ではない。

独立SOURCE295b/a5ffで全coherent q7式を既に再構成している。
今回、その式の全leafを独立root結果と既存入力から供給し、
全2673個の積項と両headの全八出力を、逆順の有理数算術で照合した。
canonical変換 `rg=(5/3)A`, `k2=B−(121/48)A` も確認した。
新しく確定した四パラメータは次の通りである。

| head | p（rg） | b（k2） |
|---|---:|---:|
| 29 | 13621/40344 | −1573/1920 |
| 30 | −32129/161376 | 336618121/20528640 |

safeなのはdeltaの差全体であり、二つのeven root個別のB63像を真の関数値とは主張しない。
zero skeletonは未知socle入力が零という仮定ではなく、全105勾配零のtheoremによりdeltaだけに適用する。
bB等の出力はplane10に限ってtrueである。
四つのsource-basis parameterが確定しても、他headのhighs、全Kbar作用・extension、
元H0やnil三乗についての結論は加えない。

独立signed replayの事前全所有上界は90,521で、明示追跡最大は72,459だった。
working polynomial・sorting key・geometric powersの一時領域は事前予備へ含めた。
全係数のbefore-operation guardを維持し、実測stored62 bit／pre112 bit、
実行0.553秒、2 GB・55秒・fresh RAW2 MB内で終了した。
失敗時のprefix保存と自動再実行禁止を備えるが、今回のsigned検証は一回で完了した。

再現コードは `code/ambient_B18_C_connecting_heads2930_values_review_20261105_v63.py`。
primary、source gate、以前の検証reportのbyteは変更していない。
