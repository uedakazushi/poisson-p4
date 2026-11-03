# 既知51個の階数四成分の S61 座標カタログ

2026年11月3日。MPS 40成分と追加11成分の一般族を、同じ61自由座標へ移した。
これは既知成分の一般点の座標化であり、全 Poisson 空間の成分上界や一般ファイバーの尽くしではない。
全57成分という下界を変えず、定義素イデアルの完全な生成系も別課題として残す。
T01 の既存の全域飽和素イデアルと T39/T40 の二次 graph 素イデアルは完了済みのまま維持する。
凍結 index の各 `defining_prime_ideal_constructed=false` はこの三件について不正確であり、
`code/data/known51_coordinate_metadata_clarification_20261103_v62.json` が現在状態を訂正する。
既存の座標・入力 hash を破壊せず、この小さい metadata 訂正を合わせて読む。

## 一般族と成分の同定

MPS 40件の構成・全パラメータ Jacobi・次元・局所 Torelli 依存は
`MPS_S61_CATALOGUE_20261103_v62.md` に記した。
追加11件の元の式は `code/additional_rank4_families_20261103_v62.py` に統一した。
その各一般係数を残した projective Jacobi を厳密に検査した。
既存の成分性の証明を入力に、群作用・全体スカラー・実際の族の微分を全150係数で作り、
Euler 鉛直24方向を除いた行列の非零小行列を新たに保存した。
有理数の分母が32003の単元である点での非零値は、標数零での階数下界だけを証明する。
上界と成分同定は表に示す既存の証明を用いる。

| ID | 射影次元 | S61次元 | 元の証明 |
|---|---:|---:|---|
| R4_DOUBLE_SMOOTH_QUADRIC |22|13|v03 `p4_poisson_completeness_v3.tex`、明示的第47成分|
| R4_TRIPLE_HYPERPLANE_QUADRIC |22|13|v05 `triple_hyperplane.tex`|
| R4_CA, R4_CB |28|19|v06 `ambient_families.tex`、Theorem AB|
| R4_H_QUARTIC_ESSENTIAL5 |23|14|v06 `nilpotent22.tex`、quarticcomponent|
| R4_ZERO_LOG_5H |23|14|v06 `zero_log_component.tex`、Theorem zero_log|
| R4_MINIMAL_M |23|14|v07 `minimal_thickness.tex`、Theorem minimal-full|
| R4_A0 |23|14|v07 `nil22_v7.tex`、A0 の全 P4 成分|
| R4_DPM |23|14|v07 `dplus_fullp4.tex`、Theorem dplus-full|
| R4_CPLUS_EXT |24|15|v08 `cplus_component.tex`、Theorem cplus-family|
| R4_CMINUS_MAG |23|14|v08 `cminus_component.tex`、Theorem cminusmag|

DPM は正確な係数体 `QQ[r]/(r²−3)` 上で扱う。
Jacobi の分子を r²−3 で割った余りが零であることを確認した。
新しい階数小行列は明記した residue embedding r mod32003 で評価している。
その後の座標写像の特殊化検査は有限体ではなく QQ(√3) の厳密演算で行った。
D+ と D− の二つの埋め込みから新成分を二つ数えてはいない。

## 全 S61 成分を優越する有理演算写像

一般族 P の上に、非零スカラー λ と一般点 p=(1:p1:p2:p3:p4) を残す。
下三角行列 L の第0列を (1,1−p1,...,1−p4) として `Q=L_*(λP)` とする。
すると Q の点 (1:1:1:1:1) における値は、元の P の p での値である。
既存の固定一点の値正規化を Q に適用し、さらに **全 Sp4 の10変数大セル** を作用させ、
Euler 商、有限 perspective 正規化、61自由係数の読み取りを行う。
反転するものは値の pivot、値の Pfaffian、Sp4 の2行列の行列式および λ である。
これらの条件は各 archive の演算節点として明記される。

この写像が全 S61 成分を優越する理由は次の通りである。
一般の PGL5 座標変換の下で、固定点 e0 の逆像 p は上のアフィン chart に入る。
p を e0 に移すことを L が与え、そこで値を J にする枠の自由度は Sp4 で尽くされる。
射影枠の共形スカラーは λ を残したことで吸収でき、残る perspective 四方向は一意に正規化される。
従って元の族の PGL 飽和が成分に稠密であるという旧定理から、本座標写像の稠密性が従う。
MPS に固有の対角トーラスを追加11件に仮定してはいない。
追加11件の源は冗長であり、源と像の次元が等しいことや一般有限性は主張しない。

出力が全パラメータで元134式を満たす証明は、一般族の Jacobi 恒等式、
有限射影変換による保存、既存の全 Poisson 模型への S61 還元定理である。
全11件で、別の旧 normalizer と全61座標を正確な特殊化で比較し、
元134式・全150復元係数・projective Jacobi を検査した。
これは実装・符号規約の確認であり、特殊化から一般恒等式を推測したものではない。

演算 archive は有理定数・変数・和・積・逆元からなる有限 DAG である。
巨大な分子を展開せずに実行できるため、有限式の構成自体を未実行の elimination に先送りしていない。
一方、これを定義素イデアルの展開生成系の完成と呼ぶこともない。

## 保存場所と再現

- 追加11族の Jacobi と新しい像小行列：`reports/ADDITIONAL_RANK4_FAMILIES_20261103_v62/`。
- 全11座標写像：`code/additional_rank4_s61_circuits_20261103_v62.py` と `reports/ADDITIONAL_RANK4_S61_CIRCUITS_20261103_v62/`。
- 全51件の index：`code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json`。
- index 生成器：`code/known_rank4_catalogue_index_20261103_v62.py`。

各 index は入力ファイル、archive、復号後の内容の SHA256 を持つ。
同じ生成器を別の新しい出力ディレクトリに実行して比較できる。
実行時間等の測定値は再実行ごとに変わるため、係数・演算・検査内容を比較する。

全定義素イデアル、特殊一次ジェットでの交叉と非被約構造、および既知51件以外の
階数四成分を除く根基比較は引き続き未解決である。

## 一次ジェットの像への接続

全51 circuit から16一次出力の微分を自動微分し、保存した正確な源の点を32003へ写して
非零小行列を作った。`reports/KNOWN51_RANK4_FIRSTJET_DIFFERENTIALS_20261103_v62/`
に全微分行列、源の列順、元の archive の hash を保存した。
T01,T39,T40 は階数16である。他の37個の MPS 成分と追加11件中 CA/CB 以外は
S61 成分自体の次元が15以下であり、16次元の一般一次ジェットを優越できない。
CA/CB の一次像は旧 `AMBIENT_WEIGHT1_IMAGE_20261030` の閉 line incidence により正確に9次元であり、
今回の微分小行列も階数9を与える。従って、既知51件の中で一般一次ジェットを優越するものは
正確に T01,T39,T40 だけである。

これに T01 の既知の次数3と二つの全 graph を合わせると、一般一次ジェット上では
**既知51成分に属する点は既知の五点に尽きる**。
これは未知成分・未知枝からの点を除いたという主張ではない。
この最後の根基比較が全一般ファイバーの尽くしの未解決部分として残る。
微分が保存点で小さい例では標本値を像次元の上界に採用せず、
別の幾何的上界と下界が一致するものだけを正確な次元と記した。
