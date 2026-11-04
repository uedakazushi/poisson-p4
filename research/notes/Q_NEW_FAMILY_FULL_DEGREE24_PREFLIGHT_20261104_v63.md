# 全24次空間の有限な基底構成案：事前検査

2026-11-04。degree24の41部分空間が単射だったことを受け、**同じ24次の全67次元空間**を有限に扱う案を整理する。高次数への自動探索、新しい NEW25 点、点評価行列は使わない。この packet は source 出力の展開や基底 rank の計算をまだ行わない。

## 一つだけの形式的 ambient 二平面

凍結済み `DEGREE16_INVARIANT_PROBE` の最初の三つの ambient β vector を使い、

\[
 \beta=\beta_0+s\beta_1+t\beta_2
\]

を固定素数32003上で形式的に扱う。これは actual NEW25 の点ではない。比較する行は s,t の**多項式係数**であり、12,16,24次でそれぞれ91,153,325行である。

全14個の C/D trace6 necklace と全36個の trace8 necklaceを候補にする。旧θ/H trace基底と同じ空間を張るとは仮定しない。次の rank gate を同じ一つの形式平面で確認できなければ止める。

1. F₄Inv₈の5列に既存N₁₂を優先して加え、trace6でrank9。
2. F₄Inv₁₂とSym²Inv₈のproduct24列がrank19。
3. trace8を加えてrank20。
4. 得た4個のprimitive12と1個のprimitive16で作る70個の24次monomialがrank67。

全候補は普遍的なQQ不変式で、denominatorは32003で単元。rank67ならcharacter上界67と一致し、選択した67列は全Inv₂₄のQQ基底となる。rank不足なら別平面、候補増殖、352個のtrace12などへ進まない。

70個の内訳は、低次数spectral35、primitive12×weight12 spectral の20、Sym²primitive12の10、primitive16×Inv₈の5である。

## 既存N₁₂の再利用と source 分母の厳密な上界

N₁₂ は F₄Inv₈を法として trace(CD²CD²) に等しく、既存rank6によりそのprimitive classは非零。これを最初のprimitive12として固定できれば、新しい source 出力は trace6三本とtrace8一本だけでよい。形式平面上でもN₁₂が必要なrankを増やすことは実際に確認する。

六変数 source の nonzero β 座標を調べると、B⁻¹を持つのは第二Cartan weightが正の座標だけで、そのweightは1以上。B⁰の座標は全て第二weight −1である。従ってweight零の次数d不変式では、B⁻¹座標の個数は高々d/2。**分母B指数は必ずd/2以下**となる。この証明はsource traceを展開せずに成立する。

既存のβ整数numeratorのsupportとℓ¹ coefficient normから、新しいsource traceの上界を得た。

|新出力|保証されるB取消し|項数上界|係数bit上界|compact JSON bytes上界|
|---|---:|---:|---:|---:|
|trace6・β次数12|6|27418|470|4990076|
|trace8・β次数16|8|80125|626|18348625|

三本＋一本の全上界33,318,853 bytesは32MiB=33,554,432以内。四本＋一本ならこの上界では収まらないため、既存N₁₂の再利用は保存予算にも重要である。選択されたwordでさらに小さな上界を計算できる。

全24次source列の共通分母は引き続きB¹²H²⁴でよい。しかし各列の理論的support上界374395は予定100k項capより大きい。従来41列の実測peak45994は参考であり、将来67列の完走保証には使わない。全source係数検査には別途、実際の選択word・出力項数を確認した上で明示的な bounded budget が必要である。

形式ambient基底jobは325×70以下の小行列で、60秒・2GB・32MiB上限を予定する。source展開はその成功後に別の許可を得た範囲でだけ行う。本packet自体は新しいrank・source出力を計算しない。

保存：`code/q_new_family_full_degree24_preflight_20261104_v63.py`、
`reports/Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63/report.json`。

入力三vector、全候補word、source weight split、character gate、全support/係数上界とhashを保持する。実行約0.25秒。
