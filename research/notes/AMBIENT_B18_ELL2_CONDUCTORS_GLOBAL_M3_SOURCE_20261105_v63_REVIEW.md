# gA・gB の全混合 m³ 包含：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`fd0d61dabb2751bfe54eb858aa3ba3774e35da46f1b3fb36e0b461a4dc07b239`。
primary `107e8c6b…`、whole parsed所有量補完V2、元の単位IFT・座標変換のsourceを固定した。
新しい関数応答・physical jet・quotient値は計算していない。

CAとCBのそれぞれの own F9座標におけるunit57解は、全次数のformal関数
f63(Y0,…,Y8)を定める。元の九変数二次cacheから、最後の二座標だけを残した
定数・一次二個・ordinary二次三個の全六係数を独立に確認し、全て零である。
従って psi(u,v)=f63(0,…,0,u,v) は `(u,v)^3` に属する。
source2のmonomialは指数vectorではなくslot multisetとして読んだ。
両unit57の固定Jacobianと保存逆行列の両積も確認した。新しいIFT解や逆行列は求めていない。

literal sourceの同定は次の通りである。

- gAはnode505938であり、CAのprefix0 unit57解で元f63 node7790の最初の七座標を零にしたもの。
- gBはnode299650であり、CBのretained indices `[7,8]` のpointed unit57解で元f63 node28539を評価したもの。

CAの受理済み独立AST証明273577は、全57式とclone関数の同一性まで確認している。
今回も元C0のprefix0・telescoping record・basis末尾505938を照合した。
CBも受理済みfull-relative peerと、literal57 bindings・57式・unit determinant・
retained `[7,8]`・root299650を照合した。
7790/28539という**制限前のf63そのもの**にorder3を割り当てる証明ではない。

両 own F9からcommon physical R0へのreadbackとその逆は、全九定数が零である。
両方のphysical行は同じ `[27,108,109,112,113,117,118,119,61]`。
実readbackの一次係数から九次Jacobianを再構成し、保存逆行列の両積と、
det CA=`3/128`、det CB=`96` を独立に確認した。
したがって両pullbackは、同じ原点における全次数の局所formal同型である。
有限二次tableを任意の全次数mapへ延長したのではなく、元のunit57で定義されたmapを使う。

局所準同型は `(u,v)^3` を `m_R0^3` へ送る。
よって **gA505938、gB299650は全physical R0で m³に属する**。
zi*r、zi*ellを含む全ordinary二次係数の消滅が従い、pureplaneだけの主張ではない。
`ell=5(s−r)` はpointedな可逆一次変換なので、この包含を保つ。ell²商へのbase changeでも保たれる。
zero7制限を元座標へのpullbackの**前**に取る点が本質である。

旧18882はcuratedな小集合だけの所有数であり、whole parse peakとは扱わない。
補完V2はsource2全体29,587、pointed9全体11,618、CB57全体2,835、
全小report・二重出力・codec・metadata余白を課して55,138<60,000とする。
独立reviewerも全parsed親objectとaliasを同時課金し、事前55,301に固定matrix検査の
明示的な小領域を加えて60,000内であることを確認した。
実before-operation最大は13 bit、実行5.256秒、55秒・2 GB・RAW32 MiB内だった。

他のg/h、nonunit quotient、opaque frontier、upper/HF・長さの値へは結論を広げない。
再現コードは `code/ambient_B18_ell2_conductors_global_m3_source_review_20261105_v63.py`。
