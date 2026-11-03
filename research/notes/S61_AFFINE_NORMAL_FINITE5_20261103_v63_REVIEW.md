# 全S61のnormal42基底上有限5生成表示：独立監査

`S61_AFFINE_NORMAL_FINITE5_20261103_v63` を、明示したbeta-only非空開上の **normal42多項式基底に関する有限5生成**、全未知残余のnormal射影の閉埋込み、および次数4以下の2796式のorigin saturationによる正確な記述として受理する。独立再検算は5.22秒で終了した。一般beta上の未知残余の空性やbeta射の有限性は受理範囲に含めない。

## 1. 元の全134式への独立な接続

主producerやその`load_fiber`をimportせず、`legacy/v02/data/symplectic_slice_reduced_equations.json` の **元61変数の全単項式** から再構成した。一次ジェットを保存T01点で固定し、45個のweight2変数を

\[
 u_{q_j}=u_{q_j}^{T01}+\sum_{k=0}^2E_{jk}z_k+
 \begin{cases}0&j<3,\\ n_{j-3}&j\ge3\end{cases}
\]

に直接置換した。演算は主側のFraction sparse実装と異なるSymPyのQQ・45変数polynomial ringで行った。

既存のplane記録は **未中心化のaffine座標** であり、定数項が零であるとは仮定しない。全45成分について保存T01点がそのaffine planeに乗ることを確認し、その平行移動後の線形部分がEに一致することを確認した。Eの最初の3行は単位行列で、残る42座標には各n_jが係数1で入るため、全45座標変換の行列式は1である。点、E、五行subtractionの係数分母は全て32003単元であった。

この変換後の134式をF_iとする。最初の五行を保持し、後の各行から保存された五行の定数係数spanを引いた **全129等式**

\[
 F_i-\sum_{l=0}^4b_{il}F_l=\sum_{j=0}^{41}A_{ij}(z,n)n_j
\]

を全係数で独立に検査した。この行変換はblock三角で対角1であり、元idealを失わない。

さらに各実際の左辺G_iから、z_kn_j係数Cとpure-normal係数Dを直接取り出した。Aの既存factorizationから再び同じ抽出をしただけではない。得られた129×126および129×945係数は、公開`cross-input.txt.gz`の全内容と **byte単位で一致** した。全G_iにpure-z項・定数項が無く、D_iの次数が1または2であることも確認した。

一般betaへの拡張には、先行の普遍三平面・五次元plane restriction・global polarization splitterの定理を引き続き依存とする。今回の再構成は凍結QQ中心を元134式に結び付ける検算であり、一般betaの巨大な全係数を新たに展開したものではない。追加の逆行列はbetaにのみ依存し、中心のminorがp単元なので同じbeta-only開へ広がる。

## 2. 126逆行列と全119070係数を整数演算で検査

主側の`fmpq_mat_inv`や`fmpq_mat_mul`を再実行せず、保存された逆行列Vとnormal係数Pを読み、独立に全分母を消した。

\[
 Z_C=\delta C_{sel},\quad Z_D=\delta D_{sel},\quad
 Z_V=\lambda V,\quad Z_P=\mu P
\]

を整数行列として作り、`fmpz_mat`で

\[
 Z_VZ_C=Z_CZ_V=\delta\lambda I_{126},\qquad
 Z_CZ_P+\mu Z_D=0
\]

を検査した。後者は **全126×945=119070係数** の等式である。入力、V、Pの全有理分母が32003単元であることも確認した。従って全126個の実際の積を

\[
 z_kn_j=P_{kj}(n),\qquad P_{kj}\in(n),\quad\deg P_{kj}\le2
\]

へ書き換えるcofactorは正確である。残る三行は `D_other+C_other*P` であり、全129行は126個のcross rewriteと三本のpure-normal次数≤2式へ可逆な定数行変換で一致する。

元の107MB出力はhashを照合して読み、公開物に再複製していない。独立整数verifierのsource、全検査log、分母bit数を含む小reportを保存した。これはmodpの単なるrankからの推測ではない。

## 3. Tiny plane blocksと実際の有限5生成

元134式の最初の五行からplane最高二次部分を再構成し、5×6の全係数が主certificateと一致することを確認した。その実際のz_k倍をpolynomial ringで作り、15×10の最高三次係数も全て照合した。微分ではなく単項式の掛け算なので、重複変数にbinomial係数を掛けていない。

選択した5×5と10×10行列について保存されたQQ逆行列の **両側の積** が単位行列であり、分母はp単元であった。よって二次最高部はrank5、三次最高部はrank10である。残す一つの二次単項式はz2²となる。

有限生成の議論にはcompletionやborder basisの平坦性を必要としない。cross rewriteはzの個数を一つ減らし、P(0)=0だからnormalを含む単項式はpure-normalへ減らせる。五本F_lもこのrewriteで `f_l(z)+K_l(n)` に変わり、K_l(0)=0である。三次の最高部全射により、純z三次を純z二次以下とnormalを含む項へ戻せる。後者はcross rewriteでpure-normalへ戻る。二次rank5でさらに五個の二次単項式を消せる。z次数による帰納により、実際の全Aは

\[
 R\{1,z_0,z_1,z_2,z_2^2\},\qquad R=K[n_0,\ldots,n_{41}]
\]

で生成される。元ideal内の正確な等式を用いたR-module有限生成であり、R上の自由性やbeta基底上の有限性ではない。

## 4. 残余閉埋込みの議論

先行のscheme-theoretic splitter `A=E5×B` と `(n)A=(1-e)A` を保持する。n=0の全fiberは実際のE5なのでB/(n)B=0である。AのR有限性からBも有限であり、原点でのNakayamaによりB_m=0となる。betaを係数基底として残す場合も、origin section上の各素idealで同じ有限加群の議論を行える。

各D(n_i)では `z_k=P_ki/n_i` からR_(n_i)→B_(n_i)が全射である。有限R-module `coker(R→B)` は、originでもその補集合でも零となるため大域的に零である。よってR→Bは全射で、SpecBはnormal空間内の閉部分schemeとなり、originを避ける。有限性を先に証明してからNakayamaと局所化のgluingを適用しており、有限体の空性を持ち上げた議論ではない。

## 5. 2796式とsaturationのexact scope

三本のH、`n_i P_kj-n_j P_ki` の3·C(42,2)=2583本、および五つのF_lを42個のchartで分母消去した210本は、それぞれ次数≤2、≤3、≤4である。

D(n_i)では互換式が全cross rewriteを戻し、210本のうち該当する五式が全F_lを戻す。従ってJの各punctured chartは **実際の全B** とscheme同型である。既にBは閉でoriginを避けるので、そのidealはm-saturatedであり、

\[
 I_B=J:(n_0,\ldots,n_{41})^\infty
\]

となる。これは根基だけの一致ではなく、未知の非被約構造も保持する。

projective incidenceについても、B上ではnormal idealが単位idealなので[v]=[n]が定まる。E5側では選択normal42行が可逆という先行証明によりprojective kernel fiberが空である。このpropernessはA45への射についてであり、betaへのpropernessではない。

Jのsaturationが単位idealになること、Bの空性・次元、betaへの有限性、一般fiber長さ5、全成分の尽くしは今回の結論ではない。これらを示すには、例えばJにnormal idealのある冪が含まれるという追加の実ideal証明が必要である。

## 保存物

独立codeは `code/s61_affine_normal_finite5_review_20261103_v63.py` と `code/s61_affine_normal_finite5_integer_review_20261103_v63.cpp`。report、整数演算logと小certificateは `reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW/` にある。主packetと先行residual42の全input/evidence hashを現行ファイルと照合した。新しいGBや一般逆行列の再計算は行っていない。
