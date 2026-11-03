# 新25次元成分の一般非被約性：三変数関数体 rank98 の独立監査

2026-11-03。`Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63` の全関数体 rank98を承認する。
実25次元成分の次元証明と合わせ、**全 Poisson scheme はこの成分に沿って一般に非被約**であり、
generic な横断 Artin 局所環の埋め込み次元は2、長さは少なくとも3となる。
正確な一般長さ5や一般二次 obstruction rank1は、この rank 証明だけからは結論しない。

producer の directional linearization を使用せず、独立の物理 Jacobi 全350二次係数辞書を作った。
各係数を原三変数 canonical150族で微分し、凍結された全150×126 canonical basisを掛けた。
QQ[n0,n1,n2] 上の全350×126係数が保存行列と一致した。
その後、70個の有理定数単元 pivot を疎多項式の別実装で再実行した。
各原 row/column番号、pivot、残る全係数を確認し、全73零行と全10零列も恒等的な零として検査した。
残る207×46の全係数は保存残余と一致する。

実行完了済み native 入力の**全**有理函数を別の parser で読み、分子・分母を一つずつ上の独立残余へ対応させた。
native log は28 pivotと成功時の全零行 assertion を伴う exact field rank28である。
既存の凍結 FLINT rational-function engine の source も点検した。
非零 field pivot の行消去だけを用い、各演算は関数体上可逆で、最後の全残余行の全 entry が零であることを assert する。
したがって70+28=98は上界・下界とも exact な QQ(n0,n1,n2) の結論である。
192.52秒の native 計算自体を今回さらに再実行したとはしない。
独立再構成と native 入力の全 entry binding、保存された実成功 receipt、既承認の exact engine を使用する。
別に指定98次小行列式を QQ の n111 点で再計算し、非零性と保存値も確認した。
この標本は下界の検算であり、上界の根拠ではない。

一般化の幾何的範囲も確認した。実三変数族の PGL・scalar 飽和は新成分に稠密であり、
canonical Jacobi rank はその作用で不変なので、三変数関数体の rank98は成分の一般点に適用される。
一般 affine tangent dimensionは126−98=28、射影では27である。
一方、成分自体の次元は独立に affine26・射影25と確定している。
他の相異なる成分との交叉は一般点を含まない。標数零の被約成分は稠密開で滑らかなので、
全 scheme がこの成分に沿って一般被約なら tangent dimensionは25になるはずで、27と矛盾する。
generic の横断 zero-dimensional 局所環に二つの余分な cotangent 方向があるため、長さは少なくとも3である。

この結論は generic 非被約性であり、被約成分自体が一般に特異であるという意味ではない。
一点 n111 の実 Artin fiber 長さ5から平坦性を仮定した議論でもない。
一般長さを5に確定するには、さらに generic 二次 obstruction の rank 上界などの exact 証明が必要である。
また単なる diagonal の n1=1 特殊化や、未検査の二変数 normal form は今回使用していない。

独立コードは `code/q_new_family_generic_jacobian_review_20261103_v63.py`、receipt は
`reports/Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63_REVIEW.json`。
