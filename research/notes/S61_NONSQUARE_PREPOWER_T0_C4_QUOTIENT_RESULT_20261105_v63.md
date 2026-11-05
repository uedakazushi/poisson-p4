# T0 の原共通分母 C⁴ 商：一次結果と独立長除算 SOURCE

2026-11-05。固定した一回の worker が完了した。
一次報告 `S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63`
の SHA は `a86f87ecddf58ae08110daf099b3fa61a0d290a018423550c22c9d68cd097131`。

\[
T_0=a^4C^4Q_0,\qquad Q_0=J_0/D_0,
\]

という正規化で、\(J_0\) は4500項、元の \(D_0\) は197 bit のままである。
全係数の \(C^4J_0=I_0\) と保存後の全係数読戻しが通った。
商の最大係数は293 bit、整数 L1 は298 bit。
実行は約0.124秒、peak53942/79558 slots、polynomial427095・scalar/root74002 products、
実入力を含む RAW は1924062 bytesであった。

この値計算は ROOT source/code review
`S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ROOT_SOURCE_REVIEW_20261105_v63`
（SHA `90be80a59a8efb499b25ca8b42f48d75a227be1582a8e2e2935e50340ed2a435`）
の受理後に一度だけ行った。準備時の誤った入力昇順 assertion は、値計算前に停止し、
`WORKER_CUSTODY..._ATTEMPT1` に code・partial input・trace・receipt を保持した。
既存7208 body の順序は変更していない。

## 別原理の独立検証 SOURCE

独立検証では7208個の元整数を remainder に直接読み、固定した \(C^4\) の25項を
四つの因子の別 multinomial 列挙から作る。a に関する先頭係数 \(100^4\) で整数長除算し、
商を求め終えるまで一次計算の4500係数を参照しない。最後に全4500係数と元分母を照合し、
remainder 全零を要求する。有限和・binomial/R-power engine の反復ではない。

入力 \(I_0\) の L1 bit 上限326と既知の Gauss 可除性から、真の商の L1 bit 上限303、
全途中 remainder 上限333を SOURCE で得る。商支持9690、remainder支持14592、固定除数25項。
expected polynomial を重複保持せず、ソート用 alias まで数えた所有量は85994、
提案 guard100000、stored518 bit、RAW8 MB、60秒・2 GB。
この SOURCE program 自体は長除算を行っていない。

\(a,C\) は FULLU 上の既知 unit なので、既存の十分条件
\(D(D_{12})\cap(D(T_0)\cup D(T_1))\) の意味は変わらない。
GCD、他の因子・tail、radical coprimality、strict cover、補集合や全 FULLU の完了は主張しない。
