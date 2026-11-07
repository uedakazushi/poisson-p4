# 単純超平面＋既約四次：純冪零十源の独立再構成

2026-11-07。課題 `E_reduced_nonnormal`。
有限 cofactor 入口の純冪零十代表について、全90係数の線形源と全残 Jacobi・実 Pfaffian を
独立に再構成した。残 Jacobi 支持の分解や既存成分への帰属は、この確認には含まれない。

入口は `HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md` とその独立確認を入力とする。
標的を単純 h と既約四次 F の積 `hF` とし、h-Casimir lift の affine 二次 bivector
`P=C+L+Q` を取る。固定する cofactor は `Z=c+Ax` である。
四つの非零 nilpotent Jordan 分割と、零又は非零 coker-translation の最大 block 長で
全十代表になる。異なる block 長の translation を捨てていない。

| 型 | Jordan 分割 | translation top 長 | 全線形 rank | 自由係数数 | 無条件の Pf 因子 |
|---:|---|---:|---:|---:|---|
| 1 | (4) | 0 | 77 | 13 | 未主張 |
| 2 | (4) | 4 | 80 | 10 | 未主張 |
| 3 | (3,1) | 0 | 71 | 19 | z |
| 4 | (3,1) | 1 | 72 | 18 | 未主張 |
| 5 | (3,1) | 3 | 75 | 15 | 未主張 |
| 6 | (2,2) | 0 | 63 | 27 | 未主張 |
| 7 | (2,2) | 2 | 70 | 20 | 未主張 |
| 8 | (2,1,1) | 0 | 54 | 36 | 未主張 |
| 9 | (2,1,1) | 1 | 62 | 28 | 未主張 |
| 10 | (2,1,1) | 2 | 65 | 25 | 未主張 |

## 算術的な独立確認

再利用した実装は `code/classification_quadratic_cofactor_review.py --nil-sources`。
producer も CAS も import せず、整数指数の簿記、標準ライブラリの Fraction 消去、
独立な AST 多項式 parser と疎多項式演算を用いた。
`jobs/hyperplane_quartic_nilpotent_sources_review.json` を gate で確認した後、
全十源の再構成が0.631秒、最大RSS 16,624 KiBで完了した。
結果は `reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_REVIEW_20261107.json`。

各型で Jordan 行列と translation を block 長から独立に組み立て、
producer の `field_linear_matrix` と `field_constant` に一致することを確認した。
全二次以下の15単項式を六 bracket に置いた90係数に対して

    div P=Z, [Z,P]=0

の全110行を直接生成した。零行も含めて保持する。
Fraction 行消去で rank と全自由列を得て、90個の解表示全てを復元し、
全110行への代入と producer の全係数との一致を確認した。
rank は標本値や部分行列の下界だけではなく、完全な線形消去の結果である。

次に、復元した六括弧そのものから全四 Jacobi 多項式を生成し、
全ての相異なる非零係数とともに producer 証明書へ照合した。
実 affine Pfaffian と、その座標次数四の最高部 `Pf(Q)` も全項再構成して一致した。
nilpotent producer の `factorized_pfaffian=null` は因数分解を実行していないという意味である。
この optional 表示だけは比較対象にせず、実 Pf と最高部 Pf は常に検証した。
既存の半単純 mode では、表示された factorized Pf との比較も引き続き行う。

型3の再構成された全 Pf 単項式は必ず z を含む。
従って全残 Jacobi を課す前から `Pf(P)=zG` である。
四次斉次化も z の因子を持ち、零となる特殊値も含め既約四次 F を与えない。
これでこの型は対象支持から無条件に排除できる。
他九型には、この検査だけで座標因子を主張していない。
因子が見つからないことを既約性の証明には使わない。

## 拡張と証明の射程

元の半単純 reviewer bytes と hash は
`reports/CLASSIFICATION_QUADRATIC_COFACTOR_REVIEWER_EXTENSION_20261107.json` に完全保存した。
拡張は nilpotent label adapter、CLI 選択、optional 因数分解表示の扱いである。
全線形行の生成、Fraction 消去、疎多項式計算と Jacobi/Pf の生成は再利用し、
旧 producer・旧証明書の bytes は変更していない。
旧 job の現行共有 reviewer binding は、この理由を明記して更新した。

確認したことは、有限入口十型の全線形源を正しく得たことと、型3の無条件排除である。
残九型の非線形 Jacobi 支持、全 irreducible quartic 点の被覆、成分数の上界、全分類は未証明。

後続の型8の全支持は `HYPERPLANE_QUARTIC_MINIMAL_NILPOTENT_SUPPORT.md` で
既約十一変数源へ尽くし、MINIMAL_Mへの全帰属を独立確認した。
本稿の全90線形源検証と、この後続の非線形支持証明を区別する。
