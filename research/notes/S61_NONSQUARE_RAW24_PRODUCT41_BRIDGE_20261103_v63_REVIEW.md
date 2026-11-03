# RAW24 から元の product41 への点別橋渡し：独立監査

2026-11-03。標数0で、元の light 変数を
`(v0,v1,v2,v3,v4,v6,v7)`、heavy 変数を `v5` とする。
元の六二次式 `q0,...,q5` は、全135式の既存の正確な weighted source reduction の
行 `(19,60,84,114,117,133)` である。今回、その既存 source certificate を依存に残し、
六式の全係数を読み直して、以下の物理座標の恒等式を独立に再構成した。
旧 producer の座標関数や行操作関数は呼んでいない。GB・追加有限体・CRT は用いない。

## 1. 元の六二次式と許された単元

元の core 座標で

```
F = 500a+90e²+3e−24,
B = (50a+18e²−9e−2)d+60e−40,
C = 100a+18e²−9e−2,   J=5ad−6e−4,
Z = 3v0−v1+v6,       W=v3.
```

開は `U={adFBCJ!=0}` である。全六式は heavy を含まず light 次数2である。
元の係数から完全除算して得る一次式 `ell` は、全係数について

```
−(de+8)q1+8q2+dq3 = W ell,
(−10a−3e²+3e)q1+3(e−1)q3+10a q4 = −a Z ell
```

を満たす。q2 と q4 をこの二行に置き換える行列式は正確に `80a` である。
従って U の各幾何点で、元の六式の張る空間は
`(q0,q1,q3,q5,W ell,Z ell)` と同じである。

ell の vi 係数を Li と置くと、独立計算で
`L6+L1=−16CJ/243` が成り立つ。よって J が非零のこの開では ell,Z,W は一次独立である。
旧 pair-cover が用いた全開の GB を再実行する必要はない。
`(v0,v2,v4,v7,Z,W,ell)` が可逆な七 light 座標であり、ell=0 上の六座標は
`(w0,...,w5)=(v0,v2,v4,v7,Z,W)` となる。
今回の実装は v1 を直接解き、v6=Z−3v0+v1 として八変数の実代入を復元した。
全座標の逆恒等式と全四二次式の制限を検査した。分母は C,J の冪と非零有理定数だけである。

## 2. Z,W を法とする独立性も全点で成り立つ

q0,q1,q3,q5 に対して W=0, v6=v1−3v0 とし、五 light 変数へ制限した。
旧 pair-cover report の選択した列を使うが、元の六式から別実装で全行を再構成し、
四つの rational pivot と全消去恒等式を照合した。pivot は順に

```
5a,   aC/240,   16C²/81,   160aC/81
```

であり、その分子・分母は既知単元だけである。従って四式は U の全点で
`(Z,W)` を法として独立である。これは関数体での generic rank だけではない。

ell 制限の公開四式に掛けた clearing multiplier も正確に
`(1,J²/25,J²/25,J/5)` である。保存された第三 core 変数は `P=ad` を意味する。
d=P/a と代入した後の a-clearing は `(0,2,2,1)` 乗であり、いずれも U 上の単元である。
`B0=(50a+18e²−9e−2)P+(60e−40)a=aB` と書けば、同じ開は
`a*P*F*B0*C*J!=0` である。新しい pair-finite packet の actual RAW24 と同じ
804係数項からの四式であり、未知因子を逆にした座標モデルではない。

## 3. 17 kernel 方向と上界41

各幾何点で light 多項式環を S、ell を法とした環を Sbar とする。
ell,Z,W が独立なので

`ell*(Z,W)*S1`

は13次元である。14本の積の唯一の関係は ZW−WZ である。
さらに四本 `ell*q0,ell*q1,ell*q3,ell*q5` は、この13次元を法として独立である。
ell を消去して `(Z,W)` を法とすれば、前節の四二次式の独立性が使えるためである。
従って元の cubic product span の ell 制限 kernel は少なくとも17次元を持つ。

RAW24 は `Sbar1*(q0,q1,q3,q5)|ell=0` の実際の24本の積である。
その階数が24なら制限像は24次元で、元の42本の積の階数は少なくとも `17+24=41` となる。
一方、全係数恒等式

```
aZ[−(de+8)q1+8q2+dq3]
+W[(−10a−3e²+3e)q1+3(e−1)q3+10a q4] = 0
```

は非零の一次 syzygy である。`W*q4` の係数が `10a` なので U のどの点でも零にはならない。
従って元の42行の階数は常に41以下であり、RAW24 rank24 の全点で正確に41となる。

同じ結論は relation を直接追っても得られる。ell を法とする24本が独立なら、六式の
一次 relation のうち四つの係数は ell の定数倍である。ell で除して Z=W=0 とすると
その四定数は零となり、残る relation は `W*(Z ell)−Z*(W ell)` だけである。
この議論は scheme の nilpotent を捨てて未知の等号を作るものではなく、各幾何点の
元の多項式係数行列の階数についての主張である。

## 4. 新しい有限不良集合定理の適用

`S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW` が、先行する gcd exception /
binary reduction の証明書と合わせ、FULLU 上の actual RAW24 rank defect を零次元以下と
独立受理したときには、この受理済み証明書を別の適用 report に hash で結び付けられる。
上の点別橋渡しから、元の六二次式×七 light 変数の product rank は FULLU の有限集合を
除く全点で41となる。元の product rank<41 の集合は RAW24 defect の部分集合である。

有限不良集合が空とは証明していない。J=0、元の full cubic rank57、quartic 被覆、
全 ambient 除去は含まない。実行済み constructive product pivot 数を増やしたとも主張しない。
**今回凍結する core report は `--bridge-only` であり、この最終適用の受理を含まない。**
新しい finite packet または全開 rank24 packet の独立監査後に、別の適用証明書を作る。
script の optional な finite 適用は独立受理 report の存在と受理 scope を確認するが、
全開 rank41をこの code が証明したと読み替えてはいけない。

再現入口は `code/s61_nonsquare_raw24_product41_bridge_review_20261103_v63.py`。
`--bridge-only` は有限不良集合の依存を使わず、点別橋渡しだけを別の work 出力へ検査する。
