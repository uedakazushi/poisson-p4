# J=0 の第二の定数単元関係：21行55列へ

2026-11-03。前の四定数単元関係で得た22行56列行列を \(M\) とする。同じ全 bivariate core 開

\[
 A=\mathbb Q[a,e,(aFB_sC(3e+2))^{-1}]
\]

で、新しい五定数単元関係を得た。\(J\) や未知の core 因子は反転していない。
\(J=5ad-6e-4=0\)、\(d=(6e+4)/(5a)\)、\(B=2B_s/(5a)\) である。

## 正確な単元恒等式

0-based の残余行13について

\[
\begin{split}
5757M_{13,7}+2400M_{13,9}+12320M_{13,26}
-22945M_{13,36}+8225M_{13,46}
=c\,a^2C(3e+2)^2,
\end{split}
\]

\[
 C=100a+18e^2-9e-2,\qquad
 c=\frac{611548345455421440945235405573910907366050899}
 {5417524312247394180878094855782598876959520000000000}\in\mathbb Q^*.
\]

native FLINT の探索が元の行に対する unit identity を厳密に確認した後、別の QQ polynomial
実装でも全係数を独立に確認した。五つの元 polynomial と差の全8項は `identity.json.gz` に
保存した。有限体や数値標本からの推測ではない。

この列和で列7を置き換える係数5757は非零有理数なので、列置換は \(A\) 上 unimodular である。
実際の witness はこの整数ベクトルの非零 QQ 倍を用い、pivot を正確に \(a^2C(3e+2)^2\) とする。
既知単元であるこの pivot を除去し、元と cokernel 同型な **21行55列**の matrix を得た。
57行の cubic 基底からの安全な unit pivot は合計36個となる。

これは残る21行の全点独立性、全余核消滅、J=0 の rank-defect support の有限性、quartic cover、
ambient nilpotence をまだ結論しない。F と \(B_s\) はこの式には現れないが、以前の消去に
使ったため、主張の基礎環では引き続き反転している。

## 実際の計算と metadata 修正

保存済み22行56列からの新しい primitive-entry span 探索で、短い7候補では指定範囲の witness
が得られず、row13で五定数 witness を得た。探索は約167秒。factor-count5 の target は古い
factor-count3 以下の探索に含まれない。

one-witness native pivot は約0.58秒で完了し、正確な21行55列を export した。その後 Python
controller の相対 path の metadata 処理で例外が発生したが、数学的な export は既に成功していた。
path を absolute に正規化する修正を行い、**探索を繰り返さず**、保存済み witness だけを別途
native replay した。出力は元の成功 export と byte 単位で完全一致した。
この fresh replay と独立 QQ 恒等式の両方を今回の report に記録した。

22行56列入力の全22行について native bivariate row gcd も厳密に調べ、全て非零 QQ 定数だった。
従って、この入力に一行だけの共通未知因子として露出する欠損曲線はない。ただし複数行の minor
による rank defect は残り得るので、全余核の消滅とは区別する。

## 保存と再現

`reports/S61_NONSQUARE_J0_SECOND_UNIT_20261103_v63/` に、恒等式、五 polynomial、witness、
unit set、探索 log、元と fresh の native pivot log、hash を保存する。元22行56列は前の
`S61_NONSQUARE_J0_PRIMITIVE_UNIT_20261103_v63` の hash に固定した。
41.1 MB の新出力は約1秒の witness replay で再生成できるため、公開保存では hash と
compact certificate を保ち、work 内の exact matrix は保持する。

小さい証明だけを fresh checkout で確認するには

```sh
python code/s61_nonsquare_J0_second_identity_replay_20261103_v63.py
```

を使う。完全な出力の再検査は
`code/s61_nonsquare_J0_second_unit_verify_20261103_v63.py` が、前の22行 matrix と保存 witness
から行う。次の探索は今回の **新21行55列**から始め、保存済みの旧 field-gcd や Smith を
繰り返さない。新しい unit が得られるたびに出力を固定して継続する。
