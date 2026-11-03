# nonsquare の J=0：正確な PID 入力と未完了の余核計算

2026-11-03。\(J=5ad-6e-4=0\)、\(aFBCd\neq0\) 上の **full cubic** 階数欠損を扱う。
本 checkpoint の確定結果は正確な入力・追加単元消去・途中証明の保存であり、余核消滅や
階数欠損集合の有限性はまだ証明していない。既存の product-rank41 欠損有限性を、full cubic
rank57 欠損有限性へ読み替えない。

## 1. J を反転しない入力

全 core 開 \(U=\{aFBC\neq0\}\) の cubic row module は自由57であり、独立検証済みの
57行基底を持つ。\(d\neq0\) 上では、原行48,55,59,57,61,22を削除する alternative ad 基底を
使用できる。続く最初の21 pivot は \(a,F,B,C\) だけの単元を使い、正確な36行70列行列を
与える。後続の \(J^{-1}\) を使う27行61列行列は今回の入力にしていない。

入力は frozen v6 handoff の `after21` であり、SHA256 は
`c65b3d80b2207e75ba794821bb41ebc3ffc864959dab27307af77cdeab9d8fd1`。
これは既存 manifest の省略大容量出力の hash と完全一致する。元57行基底の証明と21 pivot
の再現は `research/handoff/20261103_v62/v6_cubic27_coker/README.md` の最初の二段である。
今回の archive には36行70列入力そのものと21段の known-unit verifier も含めた。

\[
 F=500a+90e^2+3e-24,\quad
 B=(50a+18e^2-9e-2)d+60e-40,\quad C=100a+18e^2-9e-2.
\]

\(J=0\) では

\[
 d=\frac{6e+4}{5a},\quad
 B=\frac{2B_s}{5a},\quad
 B_s=300ae+54e^3+9e^2-24e-4.
\]

従って正確な基礎環は

\[
 A=\mathbb Q[a,e,(aFB_sC(3e+2))^{-1}].
\]

ここで \(d\neq0\) は \(3e+2\neq0\) と同値である。

## 2. 全2520成分の厳密な代入検査

native producer は各行に既知単元 \(a^{M_i}\) だけを掛け、代入後の \(a\)-poles を払う。
その行の全ての項 \(c a^p e^q d^k\) について \(M_i\geq k-p\) を選ぶ。
出力は \(\mathbb Q[a,e]\) の36行70列、108530項である。

別の Python verifier は binomial formula

\[
 a^{M_i}c a^p e^q\left(\frac{6e+4}{5a}\right)^k
 =\sum_{j=0}^k\frac{c\binom kj6^j4^{k-j}}{5^k}
 a^{p-k+M_i}e^{q+j}
\]

を用い、**全2520成分の全有理係数**を比較した。全成分一致、5.40秒。
有限体の点や数値標本から多項式恒等式を推測したものではない。

さらに既存の exact known-unit pivot engine をこの新入力に適用し、追加13 pivot を10.21秒で
得た。使用する単元は \(a,F,B_s,C,3e+2\) と非零有理数のみ。constant row/column combinations
と既知単元除算の各段を検証し、\(A\) 上で23行57列行列への正確な cokernel 同値を得た。
これは全点階数23を証明しない。fraction-free 消去によって係数量が3.8 MBから37.9 MBへ
増えたため、後の PID 計算には小さい36行70列入力を使用した。

## 3. 二つの coordinate field と、本当に必要な完了条件

\(N\) を36行70列行列の column cokernel とする。これは \(A\) 上有限生成で、その台が
full cubic の rank57 欠損集合を表す。もし

\[
 N\otimes\mathbb Q(a)=0,\qquad N\otimes\mathbb Q(e)=0
\]

を、それぞれ \(\mathbb Q(a)[e]\)、\(\mathbb Q(e)[a]\) と既知単元の localization で厳密に
示せれば、有限生成性から非零 \(h(a),k(e)\) が \(N\) を annihilate する。従って
\(\operatorname{Supp}N\subset V(h(a),k(e))\) は有限となる。

coordinate field ではその coordinate の非零多項式を係数体の単元として用いるが、もう一方の
polynomial variable の未知因子を単元として除いてはならない。片方の coordinate-field 結果
だけからも有限性は結論しない。現在は **どちらの消滅証明も得ていない**。

## 4. 実行結果と、無効な backend 出力の除外

|方法|係数環|実際の結果|
|---|---|---|
|traced Smith|\(\mathbb Q(a)[e]\)|180 CPU秒、190.19 wall秒で timeout。input-ready のみ。Smith/annihilator なし。|
|traced Smith|\(\mathbb Q(e)[a]\)|180 CPU秒、190.03 wall秒で timeout。input-ready のみ。Smith/annihilator なし。|
|localized row-gcd|\(\mathbb Q(e)[a]\)|120 CPU秒、130.02 wall秒で timeout。11段の完了は確定、12段目の witness まで保存。全余核消滅なし。|
|保存済み12段 prefix の replay/export|\(\mathbb Q(e)[a]\)|復旧後の60 CPU秒、65.01 wall秒でも final PASS/export に届かなかった。|

AS cap は各2.2 GB。timeout は数学的な unit/nonunit 証明ではない。

最初の row-gcd 実装は Singular built-in `extgcd` が rational-function denominator に対して
`conversion error: denominator!=1` を返した。poly/poly 除算でも非自明な rational-coefficient
fixture が同じ問題を再現した。この無効な実行の後半に出力された `NONUNIT_FACTOR` は、
**実際の欠損曲線を意味しない**。receipt は `arithmetic_error` に訂正し、全数学的フラグを
false とした。archive に診断 log を残すが、この候補を研究結果として使用してはならない。

現在の実装は、leading coefficient の number 比を用いる明示的な univariate long division と
extended Euclid を使う。各段で degree の減少、除算恒等式、Bezout 恒等式を検査する。
非自明な \(\mathbb Q(e)\) rational-coefficient fixture と、その保存 proof の再実行は成功した。
上表の有効な120秒実行には arithmetic error はない。

## 5. localized row-gcd の数学的な意味

active 行列 \(M\) の一行全体の gcd を \(g\) とする。未知の polynomial-variable 因子が
\(g\) に残れば、その行座標への projection は実際の cokernel から \(R/(g)\) への全射を
与えるので、欠損曲線の証明になり得る。しかし今回はそのような有効な結果を得ていない。

\(g\) が既知単元なら、Bezout column \(T\) により \((MT)_p=g\) を作る。全 pivot 行成分は
\(g\) で正確に割れるため \(H_j=M_{pj}/g\) は polynomial であり、追加 column を消去すると

\[
 M'_{ij}=M_{ij}-(MT)_iH_j\quad(i\neq p)
\]

が得られる。column をその線形結合として append することは image を変えず、pivot 行を
\(g\) で割ることは localization 上可逆である。従って各成功段は exact cokernel 同値である。
大きい70行70列 Smith transformation を保存する代わりに、各 \(T,g,p\) を保存した。

今回の有効な実行の log は12 pivot を表示し、4261 byte の proof prefix に12個の witness が
完全に保存されている。ただし log 表示と witness の書込みはその段の最終 matrix 更新より前に
行うため、**保守的には11段完了＋12段目保存済み**とする。12段全体の replay と24行70列
resume matrix の export は、この checkpoint 時点ではまだ確認できていない。

## 6. 保存と次の操作

`reports/S61_NONSQUARE_J0_PID_20261103_v63/continuation.tar.gz` は、入力、単元集合、
全係数 verifier、二つの Smith timeout、無効 backend log、有効 row-gcd
prefix、および replay/export script を元の `work/...` 相対 path で保存する。実行 binary は
含めない。37.9 MBの23行57列出力は、元 workspace に保存したまま、公開 archive には
重複収録しない。約10秒で再生成でき、`omitted_easy_to_regenerate_outputs` に hash、byte数、
正確な再生成方法を記録した。`report.json` と member hash により byte 同一性を確認できる。

新規 code は `code/s61_nonsquare_J0_pid_*_20261103_v63.{py,cpp}`。
新しい workspace では archive を repository root に展開し、既存 Cloud 設定を使う。

次は保存済み cofactor を **各段の progress marker 付きで replay** し、24行70列の正確な
matrix を export する。最初の12段の gcd 探索を繰り返す必要はない。新たな総600秒・2.2 GB
以下の budget を与えるなら、この replay と残る余核計算の continuation に割り当てる。
同じ traced Smith を変更なく再実行しない。

Cloud unavailable の間には replay 結果を読めず、復旧直後に初めて log 空、JSON・resume 行列
なしを確認した。その後の60秒 replay も未完了である。この現実の状態を、実行済み／成功済み
と仮定せず固定した。J=0 の full cubic rank-defect support の有限性、残る quartic cover と
ambient nilpotence は引き続き別の未解決問題である。
