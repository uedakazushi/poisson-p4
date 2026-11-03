# 新一次ジェット超曲面の方程式次数≥16：独立監査

2026-11-03。対象は新25次元成分の15次元primitive一次像Bであり、
全一般βファイバーの尽くしや方程式そのものは対象に含めない。
低次数baseline、degree8、degree12/14の三packetを、凍結原稿を変更せずに監査する。

## 完全上界を標本から独立に得る

閉かつprimitiveな線形二形式の全24係数について、閉性4式とprimitive4式を独自に組み立て、
16次元基底を再構成する。その三つの共変tensor添字への作用から全十表現行列を作り直す。
保存された全行列、Cartan weights、交代形式と一致する。

producerのmonomial coin-DPとは別に、Adams指標p_kとNewton恒等式

    d*h_d = sum_(k=1..d) p_k*h_(d-k)

で全整数weight multiplicitiesを次数14まで再構成する。
各次数の総数がbinomial(15+d,d)に等しいことと全除算の整数性を確認する。
さらにsigned permutationsを使うproducerとは別に、四つの正根
(2,0),(0,2),(1,1),(1,-1)のWeyl denominator積を展開して不変部分の次元を抽出する。
結果は0:1、4:1、8:5、12:9、14:3のみである。
保存された各Weyl項のweight multiplicityも全次数で照合する。
これは不変式空間の完全な上界であり、評価点の階数とは別の根拠である。

## pullback anti-representationの符号と共変性

ここでR_Aはpullback作用の微分なので、全100組について

    [R_A,R_B] = -R_[A,B]

を満たす。交代形式ΩはR_A^T Ω+ΩR_A=0を満たす。
μ_B(b)=b^T ΩR_B b/2とおくと、R_A b方向の微分は

    δ_A μ_B = μ_[A,B]

となる。trace Gramでμをsp4行列へ戻せばδ_A μ=[μ,A]であり、
theta=R_μの微分はanti符号により

    δ_A theta = [R_A,theta]

となる。従ってthetaはbの作用と同じ向きにconjugateする。
監査コードはこのmoment微分を全100組の定数行列恒等式として検証する。

F=Tr(μ²)/2は不変である。Sがb上のsymplectic作用なら
HessF(Sb)=S^(-T) HessF(b) S^(-1)なので、
H=Ω^(-1)HessFもH(Sb)=S H(b) S^(-1)を満たす。
したがってtheta/Hの全cyclic trace wordは不変である。
定数Poisson括弧grad(f)^T Ω^(-1)grad(g)も同じsymplectic作用で不変となる。
これらは標本上の不変性ではなく、全表現と上記普遍恒等式からの証明である。

## 実際の像上の下界と単射性

三producerを新しい出力先へ再実行し、全証明書をbyte比較する。
また九個の実有理点を、別のFraction回路解釈器で評価する。
回路の値を信じるだけでなく、元のcanonical斉次二ベクトルから有限射影移動・線形正規化・
S61 perspective除去を直接行い、全61出力と実primitive24係数を比較する。
全評価点でvalue pivotとPfaffianが非零である。

全九点の12次trace、最初の五点の8次traceを別に計算する。
14次非零minorに必要な最初の三点では、五個の8次不変式の全16偏微分を
前向きdual-number行列積で独立に構成し、全保存gradientとPoisson括弧に照合する。
8,12,14次の5×5、9×9、3×3有理小行列式をそれぞれ再計算する。
これらの非零性と、上で独立に得た不変空間次元5,9,3から、
各次数の全不変式のBへの制限が単射であることが従う。

低次数baselineの四次不変式についても全十Lie微分が零であること、実点で非零であることを確認する。
次数2,6,10は不変式空間自体が零である。
Sp4中心−Iは三共変添字のtensorへ−1で作用するので奇数次の不変式はない。

## 方程式次数の結論と範囲

Bは既に受理済みの既約Sp4安定な斉次超曲面である。
その主素イデアルの生成元はSp4のcharacterが自明なので不変、scalar安定性から斉次に取れる。
次数14までの非零不変式が全てBへの制限で消えないため、生成方程式の次数は少なくとも16である。
評価点は独立性の下界にのみ使い、恒等零や不変式次元の上界を標本から推論しない。
方程式の構成、次数16での存在、全高次数不変式、一般Poissonファイバーの全点は未決である。

独立コード：`code/q_new_family_invariant_degree16_review_20261103_v63.py`。
receipt：`reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json`。
三producerのfresh replayと全独立検査が成功することをreceiptの作成条件とする。
