# Q1/B の固定 S-cache worker：実行前 custody

2026-11-05。SOURCE `54b7b906…` の固定方式を実装した準備資料。
canonical4811係数とD197はそのまま複写し、headerだけを変更する。compileと`--help`のみ実行し、
Bの層・商・他因子・GCDは計算していない。

\(B=b_1d-R\)、\(b_1=50a+18e^2-9e-2\)、\(R=40-60e\)。
**b1をFULLU単元とは仮定しない。** sourceのprimitive/Gauss補題による普通のB付値を調べる。
分母はD197とsymbolic b1¹⁸。b1=0のFULLU点を捨てたり、物理的PGL変換を主張したりしない。

## 十九本の固定 cache と全 lifetime

\(S_h=R^h b_1^{18-h}\) は19本で支持合計上限3610。
まずR累乗190 slots、一つのcurrent/next b1累乗（同時324+361以下）で各Sへ直接加算する。
これらのrootを作り終えたら、**R/current/nextを全解放してから** Hasse accumulatorを確保する。
root構築34592、main37436 slots。係数とpacked keyの二セル、全metadata16000を含む。

各k=0..18について固定式
\[
P_k=\sum_{l\ge k}I_{ijl}{l\choose k}S_{l-k}a^ie^j
\]
を一accumulatorへ加算し、最初の非零で止める。層はa≤34,e≤69,total≤86−k、最大2297項。
keyはa*128+eで、全加算でe<128のためcarryはない。
逆順全項の引き戻し零と再生を同じaccumulatorで確認し、全係数を読戻す。

L1上限は301+135=436 bit（相関した二項和79^(18−l)179^lを使用）、hard518。
root構築22078とfirstnonzeroまでの全order、二repeatを含むsource上限12768337 products。
CPU完了を保証せず、45秒soft・60秒hard・2 GBで正確なprefix停止を許す。

`Q1BHP001` は全input/R/current/next/S/accumulator、D/L1、current scaledinputbase、全phase/cursor/counterを保存する。
guardはcommit前。零層ファイルも保持する。STOPは完成・非零証明ではなく自動retryを行わない。
I/O failureで二つ目の全prefixを作らない。
全state4642060 bytesを5 MBで予約し、準備入力・全history≤2.8 MBとreceipt0.2 MB込みRAW8 MB。

実際の値計算には新しいROOT source/code/whole auditとONE GOが必要である。
本資料だけでB因子、二tailの互いに素、実bad Fitting ideal、補集合消滅、全FULLU57/239は主張しない。
