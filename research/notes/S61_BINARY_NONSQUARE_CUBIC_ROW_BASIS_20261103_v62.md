# 全core開の三次row moduleは自由57：六つの厳密な関係

2026-11-03 / v62。U={aFBC≠0} の**全域**で、元135式から得た63×91の重み三matrixのrow moduleに、明示的な57行の基底を与えた。従って三次階数は全Uで57以下、一般点で57である。ここで自由なのはrow moduleであり、cokernelの自由性や全点での階数57下界を主張しない。

    F=500a+90e²+3e−24,
    B=(50a+18e²−9e−2)d+60e−40,
    C=100a+18e²−9e−2.

この結果には d≠0、J=5ad−6e−4≠0 を仮定しない。これらを使った探索表示から得た関係を元のQQ[a,e,d]へ戻し、全係数で再検査した。

## 入力と単元正規化

入力は [現在の135式から再構成した全三次matrix](S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62.md)。full入力SHA256は `dae59904925260605628d925ae0565de084c9af0ae1bccf1774f1af067830ad6` である。

63行をa,F,B,Cの共通因子と非零QQ定数だけで正規化し、正規化行をR_iと書く。検証器は元の各行がR_iと記録した単元倍率の積に等しいことを全91成分で確認する。Jやdの因子を除く必要はない。

## 六関係と二つの基底

原row番号は0始まり。保存した六関係を用い、まず五つのgenuine cubic row、次に一つのquadric product rowを順に解く。

|削除順|削除する原row|その関係での係数|
|---:|---:|---:|
|1|48|1|
|2|55|1|
|3|59|1|
|4|57|a|
|5|61|−a|
|6|14|12a/5|

各関係は、それより前に削除した行を含まない。従ってこの六列の係数matrixは三角形で、行列式は正確に −12a³/5。これはUの全点で単元である。全ての関係を正規化matrixの全91列に掛け、QQ[a,e,d]の多項式として零となることを検査する。よって六行を残りの57行で表せ、row moduleは両方向に等しい。

最後のquadric-product関係では、原row22の係数も −ad/10 である。d≠0上では、最後に14の代わりに22を削除でき、この場合の三角行列式は a³d/10。この第二の57行matrixも保存した。後続のad開では、記録した30個の単元pivot行を保持するために使える。

## 一般独立性と代数的な階数上界

第一の57行matrixを (a,e,d)=(1,2,3)、法32003で評価し、新しく選んだ57列の行列式が27632≠0となることを検査した。全有理分母がこの素数で単元であることも評価時に確認する。列選択は証明書に保存した。従ってQQ(a,e,d)でもこの57行は独立であり、R=QQ[a,e,d,(aFBC)⁻¹]上の57行の間にR係数の関係はない。生成と独立性から、これはrow moduleの自由R基底である。

全Uでの階数上界57は、今回の六つの多項式関係と単元minorから直接従う。従来の [一般長65による上界証明](S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.md) を今回の代数的証明の前提にする必要はない。

一方、自由row moduleから自由module R^91への埋込みを特殊化すると、その像の階数が落ちる可能性は残る。したがって自由row moduleを、cokernelの自由性や全特殊点での三次階数57と同一視しない。残るgenuine cubic/quartic欠損曲線の除去は引き続き別課題である。

## 再現

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_binary_nonsquare_cubic_row_basis_20261103_v62.py \
      --output /tmp/nonsquare-cubic-row-basis-new

[係数データ](../../code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json)、独立のQQ [検証器](../../code/s61_binary_nonsquare_cubic_row_basis_20261103_v62.cpp)、[producer](../../code/s61_binary_nonsquare_cubic_row_basis_20261103_v62.py)、[凍結report](../../reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json)を保存した。検証器は探索の有限体kernelを信用する必要がなく、元matrixから六関係、単元正規化、二つの三角minor、新しい一般点minorを再検査する。残る二つの57×91 matrixも圧縮して収録する。
