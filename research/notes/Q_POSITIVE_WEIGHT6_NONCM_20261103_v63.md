# 正 weight6 軌道上の全 Poisson 環は非 Cohen–Macaulay

2026-11-03。正軸の全 eta≠0 点において、**全 Poisson 局所環および完備局所環は
Cohen–Macaulayではない**。これは新25次元成分自体の被約局所環についての主張ではない。
一般長さ5やgeneric二次rankの計算も、この系には使用しない。

P+における実八変数形式横断環をAとする。
全350式と実形式sliceの独立監査により embdim A=8であり、実25次元成分の存在と上界25からdim A=3である。
完全な五次initialの独立証明は、真の全initial idealに含まれる斉次ideal Jの次元3・次数5を与える。
従ってQQ[z]/Jからgr_m Aへの斉次全射があり、両環の次元は3に等しい。
Hilbert functionを比較すれば、Hilbert–Samuel多重度は

    e(A)=e(gr_m A)≤degree(QQ[z]/J)=5

となる。低次jetの次元をactual次元の下界にした議論ではなく、別に証明されたactual次元3との比較である。

もしAがCohen–Macaulayなら、無限剰余体上で一般の三つの線形parameterを選び、
極大idealの最小reductionを作れる。このregular sequenceで割ったArtin環の長さはe(A)、
埋め込み次元は8−3=5なので、長さは少なくとも1+5=6である。
従ってCMならe(A)≥embdim A−dim A+1=6でなければならず、上のe(A)≤5と矛盾する。

23個のcanonical affine smooth因子、または22個のprojective smooth因子を戻しても、
CM性の同値とmultiplicityは変わらない。全projective環ではdim25、embdim30で同じ矛盾となる。
全非零正軸は一つのPGL軌道なので、結論はそのすべての点に及ぶ。
有限型局所環のCM性は完備化で検査できるため、未完備の局所環にも同じ結論が成り立つ。

この系は新成分の全点が非CMであるとは述べない。
n111では別の証明により全完成環がCM・type2となるため、境界の正軸上の挙動とは区別する。
被約成分の正規性・滑らかさや、ここを通る低次元成分の尽くしも今回の不等式からは結論しない。

依存receiptは `reports/Q_POSITIVE_WEIGHT6_NONCM_20261103_v63.json`、
検証コードは `code/q_positive_weight6_noncm_20261103_v63.py`。
