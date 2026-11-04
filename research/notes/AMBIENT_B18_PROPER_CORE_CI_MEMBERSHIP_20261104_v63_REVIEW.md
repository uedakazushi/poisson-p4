# 実 proper core の三生成元・CI 結論の独立監査

凍結済み PROPER_CORE_CI_MEMBERSHIP を、作者の worker・GBを実行せず、別実装の
sparse QQ 辞書および分母を払った整数係数積和だけで確認する。原 K7 は直前の独立
受理済み proper core source から読む。四つの三生成元 membership、十二の source
lift cofactor formula、全四つの元 K に対する恒等式を全係数で検証する。

unit stripping、Q0-2Q2=-Psr、最初の三生成元の ideal 変換も元 K へ結び戻す。
共通 source multiplier は M=-P D³、原点値-5であり、D=1-d² と P=d²+5 以外の
未知因子を反転しない。正の coefficient identities だけで membership は証明されるので、
保存 GB の正当性を前提にせず、GBや S-pair を再実行しない。

initial quadratic f2 の symmetric matrix determinant は144であり、f2 は幾何的に
非退化・既約である。g2 は f2 に比例しない。r=0 では g2=40a(a-v6) に対し、f2 の
二因子上の値はいずれも非零多項式。s=0 では f2=4d(-a+2d+2v6) に対し、g2 の
二因子上の値はいずれも非零多項式。この固定された四つの非消滅恒等式も別検算し、
QQ の GCD 計算だけに依存せず幾何的 coprimality を得る。

従って (f2,g2) は dimension2 の homogeneous CM 完全交叉で、r と s は共に
nonzerodivisor。initial forms (f2,g2,sr²) は次数(2,2,3)の正則列である。filtered
regular-sequence lifting と正の source membership により、実際の full K7 core は
dimension1 の完全交叉、CM・Gorenstein、multiplicity12となる。

S=QQ[[u1,...,u5,a,d,v2,v6]], J=(u5,K7) とすると、全 normalized E4 は S/dJ。
S 自体で d が regular なので 0→S/J --d→ S/dJ → S/(d)→0 が正確である。
S/J は CM dimension1、S/(d) は regular dimension8。局所コホモロジーから
H0(E4)=0、depth(E4)=1が従う。全 E4 を dimension1 や CM と主張しない。
unnormalized D0、H0(D0)、T0、full D depth にはこの結論を移さない。

実行結果・入力 hash・恒等式の検査範囲は report.json と independent-proof.json に保存する。
