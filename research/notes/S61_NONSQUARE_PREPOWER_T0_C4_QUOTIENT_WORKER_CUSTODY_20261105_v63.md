# T0 の固定 C⁴ 商：原共通分母を保つ実行前 custody

2026-11-05。これは SOURCE と worker 準備だけの記録であり、商の係数値を計算した報告ではない。

独立な ROOT 五 jet 整数 Horner 検証
`S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63`
（report SHA `60be3df6763cbfb294727209de1d2fd0bfd672b21488bd0d6753a47a83f2140f`）
により、四つの低次層が零で、第4層646係数が一次資料と一致し、正確に
\(v_C(T_0)=4\) である。元の primary と SOURCE は変更しない。

## 入力と正規化

`root239.poly` の外部 offset は4であるが、7208個の body は既に
\(U_0=T_0/a^4\) の指数を用いる。今回の入力 shift は0で、二重に4を引かない。
前の T0/C custody の整数 body をそのまま複写し、header だけを変更する。
\(U_0=I_0/D_0\)、\(D_0\) は元の197 bit の正整数である。

\[
C=100a-R(e),\qquad R(e)=2+9e-18e^2.
\]

\(C\) は primitive である（\((-1)(-9)+4(-2)=1\)）。
既知の \(C^4\) 可除性と Gauss の補題より
\(J_0=I_0/C^4\in\mathbf Z[a,e,d]\) である。
\(i=0,\ldots,16\) ごとに有限和

\[
S_i=\sum_{j\ge i+4}I_j
 {j-i-1\choose3}R^{j-i-4}100^{20-j+i}
 =100^{20}(J_0)_i
\]

を計算する構成とする。各 slab の全加算を完了する前に高い e 項を捨てない。
完成後にのみ最終支持を確認し、固定整数 \(100^{20}\) による全係数の
exact division と逆乗算を行う。内容や GCD の探索ではない。
最終保存は \(Q_0=J_0/D_0\) であり、人工的な \(100^{20}\) を分母に残さない。

最後に全係数で \(C^4J_0=I_0\) を検証する。これは有限 Laurent 展開で
残った負の a 部分が零であることも検証し、切り捨てを仮定しない。

## 実行前資源境界と停止規則

最終商の支持上限は9690、部分 slab は1387、全乗算戻しの remainder は14592。
input7208 と R cache289 に加え、係数と packed key を別々に数え、数値 metadata と
書出し scratch を16000予約する。formula phase は53148、検証 phase は79558 slots。
R cache と slab を全解放してから C の固定累乗と remainder を確保する。

scaled slab の L1 bit 上限435、固定除去後の商303、逆 scalar 乗算436、
全 C⁴ 乗算戻し333で、全整数 hard guard は518のまま。
全 scalar/root/polynomial product 上限2581132、45秒 soft、60秒 hard、2 GB。
RAW は入力・コード・receipt を含め11000000 bytes以下。全状態又は成功出力の一方だけを
保存し、I/O failure 時に二つ目の全 prefix を生成しない。

数学 guard は全 active polynomial と operation cursor を `T0C4QP001` で保存する。
STOP は商の完成や非零証明ではなく、自動 retry は行わない。
今回実行する準備 program は native の compile と `--help` のみである。
実際の商には新しい ROOT source/code audit と ONE GO が必要である。

この資料は他の tail・因子・GCD・bad locus の計算を含まず、既知の
FULLU 上の \(D(D_{12})\cap(D(T_0)\cup D(T_1))\) 十分条件を強める主張もしない。
