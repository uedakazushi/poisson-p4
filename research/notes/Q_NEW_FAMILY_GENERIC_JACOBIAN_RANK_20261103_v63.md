# 新成分の全generic Jacobi rankは98：三変数関数体の上界証明

2026年11月3日。標本(1,1,1),(1,2,3),(2,3,5)でのrank98だけからgeneric上界を推論しない。
以下は全QQ(n0,n1,n2)上の正確な行消去である。

全三変数族のcanonical150係数から、原物理Jacobiの全350 cubic係数を微分する。
独立監査済みの全150×126 canonical basisを掛けて、350×126行列Jを得る。
値・群・族の小さい行列や一部のJacobi条件で代用していない。

QQ[n0,n1,n2]における非零有理定数pivotを70回使用する。
各Schur操作とpivotを原row/column番号で保存し、全係数で再検査する。
残りの形は280×56である。そのうち**73行は全entryが恒等的に0**、
**10列も全entryが恒等的に0**なので、それらだけを削除する。
原indexと零性をすべて証明書へ保存する。残差は207×46で、全28533項である。

この全行列を既存FLINT rational-function rank engineに渡すと、
QQ(n0,n1,n2)で正確なrank28が得られた。
実際の完了出力は `EXACT_FUNCTION_FIELD_RANK 28 ROWS 207 COLUMNS 46`、
native計算192.52秒である。engineは最後の全零行の各entryを厳密にassertして終了する。
従って全generic rankは厳密に70+28=98となる。

標本(1,1,1)のQQ98×98非零minorも別に保存するが、これは下界の証拠である。
上界の証拠は全関数体の消去と全残差零行であり、標本値からの推論ではない。

generic canonical affine tangent dimensionは126−98=28、projectiveでは27である。
実際の被約成分の次元が独立に25と確定したreceiptを合わせると、この全Poisson schemeは
その成分に沿って一般に非被約である。標数零の被約既約成分は一般点で滑らかなので、
一般被約なら接空間次元25になるはずだからである。
generic transverse embedding dimensionの余剰は2である。
exact generic Artin長、二次obstructionのrank、全一般完成環は、このJacobian計算だけからは決まらない。

generic projective orbit dimension23と成分次元25の差は2であるが、
本計算では二変数のPGL normal formや有限PGL同一視を構成していない。
全三変数係数体を保った98の証明であり、n1=1などの特殊化をgeneric上界に使わない。

## 保存と再現

* `code/q_new_family_generic_jacobian_rank_20261103_v63.py`
* `reports/Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63/{report.json,certificate.json.gz,native-input.txt,native-log.txt}`

producerは全350×126を再構成し、証明書の70pivotを再検査する。
`--native` に既存rank executableを渡せばnative計算も新たに実行する。
保存時は既に実行完了した192.52秒のrunの入力と観測済みstdoutを再利用した。
その場合も再構成したnative inputとの**全byte一致**を確認し、fresh native replayとは記載しない。
保存logは成功した実際のtool stdout各chunkを連結したものである。
全数学入力・出力のSHA256を固定する。ignored workや実行バイナリを研究結果の根拠として残さない。
