# 新成分first-jet像：全16次不変式の制限単射の独立監査

2026-11-03。作者の凍結 `Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION` を独立に承認する。
新成分の既約なprimitive first-jet錐超曲面を定める斉次方程式の次数は**18以上**である。
18次の方程式の存在・非存在、制限kernel、方程式自体を求めたとは主張しない。

## 完全な不変式空間の次元

既存の独立監査が全16次元表現の三つのcovariant slot、8拘束式、symplectic形式、
全10作用行列を検査している。今回は保存された16個全ての重みを重複数ごとにまとめ、
各重みlambdaの因子 `(1-t*e^lambda)^(-m)` を二項係数 `binomial(m+k-1,k)` で展開した。
作者の個々の重みを順に加えるcoin recursionとは別の整数算法である。
各次数0--18で全重み多重度の和が `binomial(15+d,d)` と一致することを確認した。

正根 `(2,0),(0,2),(1,1),(1,-1)` のdenominator積の定数項からtrivial multiplicityを取り、
全保存Weyl項の多重度とも比較する。15,16,17,18次は正確に `(0,20,0,9)` である。
従って16次不変式空間の次元20は標本でなく、完全な表現指標からの上界・等号である。

## 全20点が実新成分の点であること

保存された18入力点を別のFraction DAG evaluatorで再評価した。
さらに元の全150係数のactual3変数Poisson族を取り、物理的projective translationと
実際の線形value normalizerを適用し、S61の全61座標がDAGと一致することを20点全てで確認した。
全value pivotとPfaffianは非零である。各primitive24座標を、保存された16次元基底に戻す等号も検査した。
従って任意に与えたW16点を実像の点と見なしたものではない。

native入力に含まれるOmega、Omega逆、Lie trace Gram逆、全10作用行列、60不変式の定義、
20選択列、20個のbeta16を一つずつ独立に再構成して照合した。
作者のQQ native処理を別出力先でfresh replayし、全4証拠ファイルのbyte一致も確認した。
独立のprime1000003上の行列算法で、20点における**全60式・全1200値**を再計算し、
保存QQ値をこの異なるprimeへ還元した値と比較した。
20次QQ行列式は別のSymPy DomainMatrix算法で正確に再計算し、nativeの長い有理数と一致し非零である。

## 不変性と次数下界の論理

thetaとHが同じ共役作用を受けることは、先行する独立reviewが全Lie作用・anti-actionの符号を含めて
普遍的に証明している。今回そのreceiptをhash固定して再使用する。
したがってtrace wordはW16全体上の真の不変式であり、各entryが二次なので長さ8は16次である。
選択された20列は全て長さ8の純粋traceである。60式の積に関する不変性は同じ普遍論証から従うが、
独立性の下界のために積列は必要ない。

非零20次minorから選択20式の実像への制限は独立である。完全空間の次元が20なので、
16次不変式の制限全体は単射となる。有限個の点での消失から恒等式を推論した論証ではない。
新像は既に独立監査された既約Sp4不変の錐超曲面であり、主素idealの既約斉次生成元は
Sp4の指標が自明であるため不変式となる。先行する14次以下の排除、15/17次空間の零性、
今回の16次制限単射を合わせると、その次数は18以上である。

一般betaの既知五点尽くし、他成分の不存在、18次のkernel、文献上の優先性は主張しない。

review code: `code/q_new_family_degree16_exclusion_review_20261103_v63.py`。
receipt: `reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63_REVIEW.json`。
