# 新三変数族の一次像は正確に15次元：一般W16を優越しない

2026年11月3日。新しい射影25次元成分の下界だけから一次ジェット優越性を推論しない。
以下では、その実際の三変数族のPGL飽和について、全一般係数体上で一次像の次元を計算する。

## 全raw jetと全source方向

affine bivector P(n0,n1,n2)を使い、評価点p=(p0,...,p3)を独立な四変数として残す。
その値C=P(p)と、平行移動後のlinear jet L=sum_i d_iP(p)*t_iを取る。
出力は六つのCと24個のL、合計30個である。

sourceの微分には三つの族方向、四つの評価点方向、GL4の全16方向、
perspectiveの全4方向、全P-scalarを含める。したがって行列は30×28である。
GL4微分は `delta C=A C+C A^T`、
`delta L=A L+L A^T−sum_i (A t)_i partial_i L`。
perspective微分はCを変えず、

    delta L=−2(alpha.t)C+t*(C alpha)^T−(C alpha)*t^T

である。これは有限projective変換t→t/(1+epsilon*alpha.t)の一次微分であり、
群方向を標本値に固定してからsource方向だけを微分する操作ではない。

## generic rankの上界と下界

QQ[n0,n1,n2,p0,...,p3]で、非零有理定数をpivotとする20回のSchur消去を行う。
全操作は多項式環上で可逆であり、全pivotと残ったrow/columnを保存する。
残る行列は10×8、全1345 monomial項である。

既存の凍結native engine `rational_function_matrix_rank_20261103_v62.cpp` へ
この全多項式行列を渡す。同engineはFLINTの有理関数体
QQ(n0,n1,n2,p0,...,p3)で可逆な行消去をし、全残り零行の各entryを厳密に検査する。
そのrankは正確に5である。従って元のgeneric raw-jet rankは20+5=25。
有限体や有限個の点の低rankを上界に使っていない。

別にn=(1,1,1),p=0で正確な25×25有理minorを計算し、非零であることを確認する。
ここで値Pf(C)=−16/81なので正規化は有効である。
このminorは下界の独立な証拠であり、上界は上の全関数体計算による。

## canonical primitive betaへの移行

非退化な値Cを標準Jへ動かす自由度は6次元であり、
perspectiveで除くtrace方向は4次元である。
同じsource行列でGL4の値部分rank6、perspective部分rank4（値部分は零）を検査する。
GL4とperspectiveの全群方向を残しているため、raw imageは局所的にこれら10 nuisance方向と
canonical primitive beta imageの積となる。
これは既存T38の全raw-jet法と同じ値正規化・perspective quotientの議論である。

従って新族のPGL飽和のprimitive一次像閉包は正確に

    25−6−4=15次元

であり、W16のhypersurfaceである。
GL飽和が新しい25次元成分全体に稠密であることを、別のactual local upper25と合わせて確定した後は、
この結論はその全成分の一次像の次元である。
S61の成分次元16に対する一般fiber次元は16−15=1となる。

この成分は一般W16を優越しない。
したがって新成分が得られたことだけから、十分一般な一次ジェットの既知五点予想や
実際の一般fiberの尽くし上げが反証されたとは言えない。
15次元の特殊一次像上の追加fiberは、一般W16の問題と区別して研究する。
像hypersurfaceの展開生成多項式や特殊fiberの全scheme型は、このrank証明とは別課題である。

## 再現と保存

既存native sourceを通常の `g++ -O2 -std=c++17 ... -lflint -lgmp` でignored workへcompileし、
`q_new_family_firstjet_image_20261103_v63.py` の `--native` へ渡す。
`--work` は一時行列、`--output` は新規report先を指定する。

* `code/q_new_family_firstjet_image_20261103_v63.py`
* `reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63/{report.json,certificate.json.gz,native-input.txt,native-log.txt}`

全30×28行列、20unit操作、10×8残差行列、nativeの全入力と実際の完了log、25minorを保存する。
