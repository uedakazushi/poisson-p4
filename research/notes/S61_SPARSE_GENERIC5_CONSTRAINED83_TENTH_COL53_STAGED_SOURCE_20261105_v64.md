# 第10段 col53 の共有SOURCEと正確な packed5 支持

全167×14の残り14列について、既存第9段の固定literal recipeの比較からcol56を除くとcol53が最小である。これは指定recipeの構造上界による選択で、実際の係数支持の全列最適性やランクを主張しない。旧raw入力・既存col82/col56値・F9/F10だけを読み、新しい行列係数は求めない。

共有N8_35,53とN9_58,53をそれぞれ一度だけ組み立てる形式DAGを用いる。U=Q^4 R^2 E、C10=353894400000 H t0として

N10_i53 = U F10 N9_i53 − C10 N9_i56 N9_58,53,
D10 = U^2 F9 F10

を維持する。全14出力以外の13列は型付きのまま、全167行・87追加RHS・66零制約RHS・P290一回・全復元式とF9/F10を含む全補集合を保持する。F9の追加除算やGCDを行わない。

SOURCEの支持だけを pack(e)=BIAS+Σe_j256^j、BIAS=64Σ256^jで正確に保存する。各指数が[-64,191]内ならbiased5座標は一意のbase256 digitであり、特殊化ではない。各Minkowski積より前に両入力BOXとshiftから出力BOXを検査し、packed和はpack(e)+pack(f)+rawindex(shift)−BIASとする。biasを引かない誤った和を使わない。支持hashはdecoded5tupleの通常のlex順JSONと一致させ、alias heapsortの参照配列一つだけを保持する。

SOURCE ownedは全source metadata・index・一時支持・hash aliasを課金し、150000で先行停止する。固定READ集合18MiB以下とSOURCE全出力8MiB以下（無圧縮）の予算を支持計算の前に保存する。60秒/2GB。各完了出力のprefixと失敗時のliteral DAG/全profileを保存し、自動再試行しない。

将来VALUEの係数写像は引き続き5指数+整数+slot等の通常のscalar数を用いて別に予算化する。SOURCE packed表現をVALUEの150k回避に使わない。共通scalar denominator・正係数majorant・全last-useとcodec RAWを計算し、予算が不成立なら係数GOを出さない。direct billion tuple展開、全行列係数、全逆行列、exhaustionは未実行である。
