# J=0 の新しい四定数単元関係：全 core 開で22行56列へ

2026-11-03。\(J=5ad-6e-4=0\)、\(aFBCd\neq0\) の full cubic 欠損問題で、新しい exact
known-unit pivot を一つ得た。\(J\) や未知の core 因子は反転していない。係数体の一方の
coordinate を一般化した結果ではなく、次の **全 bivariate core 開**で成立する。

\[
 A=\mathbb Q[a,e,(aFB_sC(3e+2))^{-1}],\quad
 B_s=300ae+54e^3+9e^2-24e-4,
\]
\[
 F=500a+90e^2+3e-24,\qquad C=100a+18e^2-9e-2.
\]

\(d=(6e+4)/(5a)\)、\(B=2B_s/(5a)\) である。

## 単元恒等式

先の `S61_NONSQUARE_J0_PID_20261103_v63` における、21個＋13個の安全な単元 pivot 後の
23行57列行列を \(M\) とする。保存済み hash によって同じ行列であることを固定する。
番号は0-basedであり、**残余行列の第7行**について

\[
 M_{7,4}-M_{7,15}-M_{7,25}+M_{7,30}
 =c\,aC(3e+2)^3,
\]

\[
 c=\frac{26807011153965784462597440300438824677423}
 {5043161936797200128766360308089875274986525000000000}\in\mathbb Q^*
\]

が成立する。四つの元 polynomial と差の全10項を `identity.json.gz` に保存した。
native FLINT の探索・元行の再検証とは独立に、Python/SymPy の QQ polynomial arithmetic で
**全係数の恒等式**を新たに確認した。数値点や有限体からの推測ではない。

右辺は \(A\) の単元である。列4を上の \((+1,-1,-1,+1)\) の列和で置き換える操作は、その
列自身の係数が1なので unimodular である。native replay では同じ四係数を \(1/c\) 倍した
witness を用い、pivot そのものを \(aC(3e+2)^3\) とした。非零 QQ 倍は依然可逆である。

従ってこの単元 pivot を除去し、元の cokernel と \(A\) 上同型な **22行56列**の正確な残余を
得た。57行基底からの安全な pivot は合計35個である。最後の22行の全点独立性、余核消滅、
J=0 の欠損集合有限性や ambient 冪零性は、この結果だけでは証明していない。

## 探索方法と検証

既存の V9 primitive-entry span 方法を再使用した。各 entry から既知の単元因子だけを
除くことは探索のために行い、その後で分母を既知単元の積により払い、**元の行**に対する
polynomial witness と単元性を厳密に再確認する。entry ごとの除算後の matrix を、元 matrix
と同値だと仮定していない。

短い4行だけを対象にした300秒以下の bounded search で、最初の3行では指定された定数 span
の witness が得られず、第4候補の row7 で得られた。今回の witness は四係数とも定数であり、
active entries の denominator clearing は不要だった。
target は既知単元因子を5個含むため、従来の factor-count3 以下の target 探索では見つからない。
この bounded refusal は、各行や全 matrix の単元 span が存在しないという証明ではない。

探索全体と exact pivot replay は65.81秒、pivot 部分は4.68秒であった。PID/Smith の
再実行、保存済み12段の field-gcd 再探索、大きい GB は用いていない。

## 証拠と再現

`reports/S61_NONSQUARE_J0_PRIMITIVE_UNIT_20261103_v63/` に、独立恒等式、四 polynomial、
witness、native search source snapshot、単元集合、search/pivot log、hash を保存する。
新しい39.5 MB残余は元 workspace に保持し、公開保存では hash と再現手順を固定する。
前の13 pivot の約10秒の replay に続き、今回の witness だけなら約5秒で再生成できる。

独立 verifier は

```sh
python code/s61_nonsquare_J0_primitive_unit_verify_20261103_v63.py \
  --output work/nonsquare-J0-primitive-identity-replay
```

である。full input は前の J0 checkpoint を展開し、保存済みの native known-unit engine で
23行57列を再生成して用いる。今回の探索を全て繰り返す必要はなく、保存済み row7 witness を
適用すればよい。新しい scalar/known-unit factor を除くときには、その exact division と
元 polynomial identity を引き続き検査する。

この新しい22行56列から、同じく bounded な方法で残る単元・rank 条件を調べることが次の
作業である。全 J=0 の cubic 欠損有限性と、その後の quartic/ambient の除去は未完了である。
