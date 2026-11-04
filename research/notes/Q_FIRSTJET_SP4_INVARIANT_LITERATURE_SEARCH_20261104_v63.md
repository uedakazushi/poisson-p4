# 16次元 Sp4 表現の不変式環：限定した一次文献検索

2026-11-04 UTC。既存の不変式計算を置き換える公式を探すための検索記録であり、
文献の不存在・優先権・新規性を主張するものではない。新しい点評価、生成元探索、
Hilbert 級数の追加計算は行っていない。

## 対象と結果

対象は C2 の最高ウェイト ω1+ω2（ε 座標で (2,1)）の16次元既約表現 W。
既受理の `Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63` により
W=ker(E5⊗S4→S4) という gamma-traceless vector-spinor と厳密に同定されている。
単なる4×4行列全体、adjoint10次元、有限群Sp4(q)、別のSpin群の16次元spinorとは区別した。

今回の検索では、**この W の全 Hilbert 級数、全最小生成元、全不変式環表示を
一次本文で確認できる文献は見つからなかった**。これは限定検索の結果だけであり、
そのような文献が存在しないという結論ではない。既存の商次元6、有限次数の
character 計算、degree≥24 の内部証明はそのまま使い、外部の生成次数上界を追加しない。

## 本文または著者原稿を確認した近接文献

1. Hanspeter Kraft and Gerald W. Schwarz,
   [Representations with a Reduced Null Cone](https://arxiv.org/html/1112.3634v1),
   arXiv:1112.3634v1。Introduction と §7 の対象を確認した。
   classical adjoint group の分類は W に直接適用できない。
   W では中心−Iが spinor 因子に−1として作用するため、adjoint group の表現ではない。
   本稿から W の coreducedness、cofreeness、Hilbert 級数を推定しない。

2. Shouhei Ma,
   [Rationality of fields of invariants for some representations of SL2×SL2](https://arxiv.org/html/1202.2947v2),
   Compositio Math. 149 (2013), 1225–1234,
   [DOI 10.1112/S0010437X13007069](https://doi.org/10.1112/S0010437X13007069)。
   §1 と §2.1 を確認した。定理は単独の V(a,b)=SymᵃA⊗SymᵇB の射影商を扱う。
   受理済み Spin4 切片は V(2,1)⊕V(1,2) という12次元直和なので、同じ表現ではない。
   一般の bundle/no-name の道具は参考になるが、almost-free などの仮定を別に
   確認せず、この直和やその normalizer 商の有理性へ移さない。

3. Peter J. Olver,
   [Invariant theory of biforms](https://www-users.cse.umn.edu/~olver/a_/biform.pdf),
   University of Minnesota preprint, August 1986。
   [著者の公開原稿一覧](https://www-users.cse.umn.edu/~olver/paper.html) から15ページの
   scan を開き、原稿の pp.1–2,9,11–15 を画像で確認した。
   Theorem 6 の反復 transvectant による生成法と、(1,1)、(2,2) の例は参考になる。
   今回読んだ範囲は直和 V(2,1)⊕V(1,2) の Hilbert 級数や生成次数表を与えていない。
   この古い原稿の存在と一般的方法だけを記録し、対象直和の既知表示としては引用しない。

## 本文の確認ができなかった文献と除外

Gerald W. Schwarz, *Representations of simple Lie groups with regular rings of invariants*,
Invent. Math.49 (1978), no.2,167–191 は Kraft–Schwarz の参考文献で書誌を確認した。
[DOI 10.1007/BF01403085](https://doi.org/10.1007/BF01403085) の本文は今回の接続では取得できなかった。
したがってその表から W の記載・非記載を主張しない。

Olver の書籍 *Classical Invariant Theory*（1999）の Cambridge Chapter10
[Multivariate Polynomials](https://doi.org/10.1017/CBO9780511623660.011) は概要のみ確認できた。
章本文の具体的な生成元は未確認で、研究結果の根拠にしない。

検索に現れた conjugation-invariant matrix algebra、Sp4 の有限群、外積代数、
symplectic quotient の Hilbert 級数は、対象表現の一致を確認できないため採用しなかった。

## 次に使える範囲

実際の12次元切片 U∈V(2,1), V∈V(1,2) と全6個の Lie 作用は既に内部証明書で
受理済みである。この正確な直和を指定した simultaneous biform transvectant の
計算や、Schwarz 原論文の表・後続論文の確認が次の文献作業となる。
一般の単独 biform の有理性や生成定理から、NEW25 の式、次数24の存在、
残余曲線の種数、全不変式環生成を導くことはしない。

主要な検索語とアクセス結果、既受理の表現同定への hash は同名 report に保存した。
文献の PDF を転載・再配布しておらず、著者への連絡や外部投稿も行っていない。
