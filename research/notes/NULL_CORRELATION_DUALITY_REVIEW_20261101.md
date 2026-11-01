# Null-correlation 双対曲線の独立監査

2026-11-01。対象は `NULL_CORRELATION_DUALITY_20261101.md` と対応コード。
結論を受理する。一般関数体上の torsor の区別と、無限遠からアフィン解を
推論しない範囲設定が必要であり、対象ノートはこれらを保持している。

## 行列恒等式の別実装

`code/null_correlation_duality_review_20261101.py` は S61 の `primitive_beta` を
使わず、24係数の closed/primitive 制約の核から16変数の β を構成した。
対象の Pfaffian helper も使わず、各添字を除いた四添字から五つの Pfaffian を
直接展開した。次のすべてが全16変数で一致した。

- d(ι_Eβ)=3β。
- Φ_E=[[β+tJ,Jy],[-(Jy)ᵀ,0]] の五つの Pfaffian。
- (Jβ)²=−Pf(β)I と、第五式が局所的に最初の四式から従う行列恒等式。
- y=0 の射影中心が Pfaffian locus に属さないこと。
- Kβ と Φ_E の全係数の直交性、および既存 Φ_C との congruence。
- 両線形パラメータ空間が常に5次元であること。

証明書 `reports/NULL_CORRELATION_DUALITY_REVIEW_20261101.json` は PASS。
同じ closed/primitive 制約の定義は共有するが、対象の高水準構成は再利用していない。

## スキーム・束・原典の確認

Jy のいずれかの係数が可逆な chart 上で t=(βy)i/(Jy)i が一意に定まる。
残る比例条件は null-correlation quotient の section の零イデアルそのものである。
第五 Pfaffian もこの局所イデアルに属するため、これは点集合の比較ではなく、
任意の β に対する零スキームと Pfaffian locus の同型を与える。
中心を含まないことも t²=0 を射影 chart t=1 で読むことでスキームとして確認できる。

束の議論では N(−2) の acyclicity を確認した。実際
0→O(−3)→Ω¹(−1)→N(−2)→0 と Euler sequence の各行束はこの範囲で acyclic。
Koszul sequence の H²(O(−4))=0 と合わせて H¹(I_E)=0 となり、
一般零曲線の連結性が従う。c₂=5H²、adjunction、大域生成性から得る
次数5・種数1・一般滑らかさの論証に欠落は見つからなかった。
W16→H⁰N(2) の単射性では、primitive trace の非零係数5/3が
quotient の4次元 kernel を正確に排除する。

原典 [Shinder–Zhang, Proposition 2.8](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/blms.12339)
の仮定・主張・証明を本文で確認した。同次元の直交する二つの Grassmannian
線形切断の一方が滑らかな横断曲線なら、他方も滑らかで Pic³ に同型となる。
既存 Cβ の Δ≠0 での滑らかな Pfaffian 曲線という結果と、上記の全係数の
直交性がこの仮定を満たす。従って Eβ≅Pic³(Cβ) と Jacobian の同型を受理する。
原点のない一般 torsor 自体を Cβ と同一視することはこの定理からは従わない。

G6 の同変性については既存
`S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.md` §6 の論証を再確認した。
syzygy 商の自明な一次元表現と Bianchi の誤差の処理は flag に制限する前の
最高次イデアルに関するものであり、135本の最高二次式全体への適用は正しい。
これは既存証明の範囲の監査であり、核36次元の別計算ではない。

今回 Fisher の全収縮は再実装していない。双対性の普遍的な根拠は
有限個の c4,c6 一致ではなく、全係数の直交性と引用定理である。
