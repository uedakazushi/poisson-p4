# 新25次元成分と既知57成分の短い内在的区別

2026年11月3日。新しい25次元irreducible saturationが実際のfull componentであることは
`Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW` の受理済み結果を使う。
本稿はそのdistinctnessの、actual divisorによる短い別証明である。

quinticの空間内で

    Sigma32={ [L³Q] : 0≠L∈V*, 0≠Q∈Sym²V* }

はprojective spacesの積からのproper imageなので閉である。
rank4 Poisson locus上のactual-quintic mapはregularである。
新成分の稠密な明示族ではD=h³Q、PGL移動後もL³Qとなるため、
**新成分の全rank4点のactual divisorはSigma32に入る**。

各MPS40 componentには、actual divisorが五つの相異なるcoordinate hyperplanesの和であるactual toric memberがある。
保存された各biresidue exampleの4×4block BについてPf(B)≠0を全40件QQで確認した。
period basis内の係数も解き、実際に同じdiagramのtoric periodであることを検査した。
A=−Pf(B)*B^−1とすると、toric tensorのactual quinticは

    Pf(A)*x0*x1*x2*x3*x4 = Pf(B)*x0*x1*x2*x3*x4

である。これはrank4で、linear factorの重複度が全て1なのでSigma32の外にある。
したがって新componentがどのMPS componentとも等しくなることはない。
MPS componentへのtoric pointの帰属と40件の相互distinctnessには、既存のMPS theorem/crosswalkを使う。

残りの11個の既知rank4成分の射影次元は22,23,24,28のいずれかであり、25のものはない。
従って新25次元componentはそれらとも異なる。
6個のgeneric rank2 full componentはrank≤2の閉集合に含まれるため、generic rank4の新componentとは異なる。

よって既存の57個と独立な新成分を得て、全Poisson空間の成分数下界は58となる。
この議論は全componentの尽くしを含まない。

* `code/r4_positive_weight6_distinctness_20261103_v63.py`
* `reports/R4_POSITIVE_WEIGHT6_DISTINCTNESS_20261103_v63/report.json`

全40toric witnessesと11次元listをQQで保存した。新規計算0.22秒。
既知11成分のactual divisorを全て5Hと仮定する議論は使わない。
