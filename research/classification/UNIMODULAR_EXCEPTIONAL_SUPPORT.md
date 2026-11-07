# E(3)最高部：主軌道は空、実境界は閉nilpotent incidenceで被覆

## 1. 主軌道の全30係数

[Pym, arXiv:1403.6444v1, §3.6, 式(13)](https://arxiv.org/html/1403.6444#S3.SS6)
の最高部を、座標順を (x,y,z,w)、bivector順を (12,13,14,23,24,34) として

\[
 Q=(5x^2,5xy,5xz,y^2+3xz,yz+7xw,7yw-3z^2)
\]

とする。これは線形場

\[
 X=(-5x\partial_x-y\partial_y+3z\partial_z+7w\partial_w)/4,
 \qquad Y=4x\partial_y+4y\partial_z+4z\partial_w
\]

に対して (Q=-X\wedge Y) であり、
\(\operatorname{div}Q=0\)、\(\operatorname{Pf}Q=0\) である。
以下では全24一次係数 (L) と全6定数係数 (C) を保持する。

**命題。** (P=Q+L+C) がunimodularで、Pfaffianが定数なら、その定数は零である。
従ってGで必要な非零定数Pfaffianを持つ延長は、E(3)主軌道の全点上で存在しない。

**証明。** \(\operatorname{div}L=0\) と全三次Pfaffian係数の行列 (M) は
24列、階数20である。平行移動の一次変化
\(a\mapsto dQ(a)\) は階数4で、\(M dQ(a)=0\)。よって全核はこの平行移動像に等しい。
ある平行移動で (L=0) にできる。元のPfaffianが定数なので、この操作はその非零性も保存する。
正規化後、全二次Pfaffian係数の定数 (C) に関する行列は6列、階数6であり、(C=0) が従う。
従って \(\operatorname{Pf}P=0\)。GL4変換は議論と定数Pfaffianの非零性を保つので、代表だけでなく主軌道全体が排除される。□

全行列、核、三つの非零minor、逆向きの平行移動補完、全残式は
[producer証明書](../../reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_20261107.json) に保存した。
[独立Fraction証明書](../../reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_REVIEW_20261107.json)
は二つの線形場からQを再構成し、全24列・全6列と核の完全性を別実装で検証している。
producerの再実行を独立検証とは呼んでいない。

## 2. 主軌道からの外挿を避ける実境界の証明

主軌道上の空性だけから、その閉包上の全延長の空性は従わない。
境界で (M) の階数が落ち、新しい低次源が出る可能性を残す。

[Alcántara–Cerveau, arXiv:2608.26470v1, §4–5](https://arxiv.org/html/2608.26470#S4)
は、葉層空間内の境界四成分と、その稠密部分の対数型三軌道・pullback族を主張する。
仮にこの全境界の主張を入力にすれば、必要な帰着は
\(\partial E(3)\subset\overline{L(1,1,2)}\cup\mathrm{PB}\)
である。実対数源のみの計算では足りず、閉L112の全支持結果が必要である。
さらに同論文の葉層空間は特異集合の余次元が二以上の開集合なので、
共通因子を持つ1-formの補集合は別に保持する。

ただし同論文§4は開軌道の抽象的積分解から、固定表現内の軌道閉包が対応する
積コンパクト化と同型であると進んでいる。開集合の同型だけでは境界での同型は従わない。
associated bundleのコンパクト化と実際の軌道閉包の間にも、一般には縮約があり得る。
この点を本研究で独立に補完していないため、**同論文の四成分記述を今回の被覆証明として採用しない**。

一方、Borel群 (B\subset G=SL_4\) に対する弱い還元
\[
 \overline{G\cdot[\theta]}=G\cdot\overline{B\cdot[\theta]}
\]
は独立に成立する。実際、(Y=\overline{B\cdot[\theta]}) はprojectiveであり、
associated bundle (G\times_B Y\) もprojectiveである。
自然写像 \([g,y]\mapsto gy\) の像は閉、主軌道を含み、軌道閉包に含まれるので等号を得る。
ここで自然写像が同型であるとは主張しない。

この固定表現内の実際の全境界を既済全源へ送る義務は、
`UNIMODULAR_EXCEPTIONAL_NILPOTENT_INCIDENCE.md` で独立に閉じた。
非零nilpotent線形場Yの射影空間で、`Y∧Q=0`、`Lie_YQ=0` を課す閉incidenceが
全主軌道とその実閉包を含む。全四Jordan型の全80三次一形式係数を処理し、
rank二の二型の完全核各七次元は稠密L112と閉性から全て閉L112へ入る。
regular型の完全核三次元は主軌道又は一次Casimir、rank一は全pullbackとなる。
全行列と完全核の算術を別のFraction実装で再構成した。
従って主軌道空性と既済の閉L112・一次Casimir・pullback全延長定理により、
E3実閉包全体の一般primitive底への非優越が証明された。
有限個の対角1PS極限や上記外部積コンパクト化は、この新証明の入力ではない。

## 3. 再現

```sh
python scripts/classification.py gate research/classification/jobs/unimodular_exceptional.json
python code/classification_unimodular_exceptional.py
python scripts/classification.py gate research/classification/jobs/unimodular_exceptional_review.json
python code/classification_unimodular_exceptional_review.py
```

対数源側の独立全16型と境界全10型の係数確認は、それぞれ
`CLASSIFICATION_UNIMODULAR_LOGARITHMIC_REVIEW_20261107.json` と
`CLASSIFICATION_UNIMODULAR_LOGARITHMIC_BOUNDARY_REVIEW_20261107.json` に保存した。
それらの支持・配置幾何は `UNIMODULAR_LOGARITHMIC_SUPPORT.md` の全枝に従う。
