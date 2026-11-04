# 実Tのcomplete5×825 Schur operator：独立全環監査

2026-11-04。先行のcomplete116×936提示、X0、W0およびlinear/二次peerを固定する。
新しいscalar coefficientやsource、rank、348行列は評価しない。

全768 sourceを覆う最初の27+741 constant column操作は、元source順で各列の対角係数が非零、
後続sourceを使わないtriangular matrixとpermutationであり可逆である。
自由32座標の五quotient rowを独立left annihilatorとして再構築し、全27定数列を殺すこと、
tail五座標にidentityとなることを確認する。27pivot rowと合わせて全32座標を保持する。

全741列への第二の84+657操作も、既受理の全column factorizationから別の三角式で再構築する。
全normalized Q1と残余linear関係の参照先を照合する。
actual27 blockの行列式−1、actual84 blockの行列式+1は、**正規化後の原点行列の値**である。
全形式関数の行列式が全次数で一定とは主張しない。原点非零だけで完成環上のunit inverseが存在し、
Schur操作にはそれで十分である。

4155 typed node、全形状・参照・constant組合せ、4125 scalar readback rootを検査する。
元sourceは実full Tの全936列を保持し、d3の168列を途中で落とさない。
一段目は `What=Wr-Wp A^-1 C`、`Yhat=Y-B A^-1 C`。
二段目は `(-Ystar U^-1 D | Yrest-Ystar U^-1 Wrest)`。
積順を保つ自由wordのunit identityでpivot blockの零化も独立検査する。
二つの可逆target/source操作とunit対消去はcokernelを全環上で保存する。
したがって全関係を持つ実5×825表示であり、関係825本の極小性は主張しない。

B,C,YとYstar,Yrest,Wrestはそれぞれ最大idealに入り、d3も定数零なので、最後の全関係は定数零。
従って五生成は極小である。`Yhat_1=Y_1`、`What_0=Wr_0`、sourceのlinearはYrest1、
d3の二次は `-Ystar1 U0^-1 D1` と、先行peerの実低次数interfaceにも一致する。
Eulerは元の形式関数を参照し、I²や有限商で先に微分しない。

T全長・K torsion・原点深さ・全nil指数は未決。
専用checker `code/ambient_B18_T_minimal5_Schur_operator_review_20261104_v63.py` と同名REVIEWに保存する。
