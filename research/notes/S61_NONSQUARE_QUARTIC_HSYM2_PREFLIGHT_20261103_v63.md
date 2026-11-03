# Genuine quartic のH*Sym²係数だけを見る bounded preflight

2026-11-03。凍結済みactual135 lazy oracleから、全33genuine quarticのH*x_i*x_j (i<=j, 7light) 係数を924個queryする。既存H²unitの2×2source変換と、他31行からのH²消去をDAG上で適用する。必要なH²係数33個もqueryする。元oracleを変更せず、新しい出力だけにcoefficient circuitを保存する。

この作業は33×28およびH²消去後32×28の**記号的support・直接row relationのpreflight**である。core parameter polynomialを展開せず、実core点・有限体・GB・rank計算を使わない。

係数nodeを独立なcore-atom記号の多項式として正規化し、零になる場合だけ厳密なQQ恒等式として認定する。atom間の実rational-function関係は無視するので、零にならなかった式を非零とは認定しない。従って保存するsupportは全て**supportの上界**であり、generic supportやmatrix rankとは呼ばない。同じformal normalizationから得た直接QQ比例rowは厳密identityである。

正規化は1node2000項、cache総250000項に制限し、拒否は未計算と記録する。DAG200000node、artifact32MiB、実行60秒・2GiBを上限とする。元oracleのhashを開始時と終了時に照合する。

H²source変換は全aFBC開で可逆である。一方、ここではまだheavy7 graphやpure146 quotientを適用していない。純四次の残余や実C122 rank64の結論は得ない。source上の疎なパターンや直接恒等式が見つかれば、次の実quartic reductionに用いる。
