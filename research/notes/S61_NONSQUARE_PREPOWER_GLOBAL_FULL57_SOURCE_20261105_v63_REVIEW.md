# FULLU 全域の old25 → current15 → actual full57 SOURCE 橋渡し：独立監査

2026-11-05。PRIMARY
`S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63/report.json`
SHA256 `f56e3d125e18ff1aa1c3d8720dc574062ff5a54759ec0bf13d17fc17267b7499`
を独立に監査した。REVIEW report は
`reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW/report.json`
SHA256 `6640c4f2be0e409a8f4854fd458fbbd700b4e28e7d6a360c153b4f2f67bcd25b`。
author producer の import、新しい行列係数・行列式・rank・点・標数の計算は行っていない。

係数環は初めから

\[
R=\mathbb Q[a,e,d,(adFBCJ)^{-1}],\quad
F=500a+90e^2+3e-24,
\]

\[
B=(50a+18e^2-9e-2)d+60e-40,\quad
C=100a+18e^2-9e-2,\quad J=5ad-6e-4
\]

である。以下の prefix では Δ₄₀₇、Saux、N6、N7bar、P_K92、h22、p、R434、
κ3、anchor A の逆数を導入していない。

## 実際の product41 と選択41行

`PRODUCT41_FULL_OPEN` を名前だけで採用せず、独立受理済み
`PAIR_EMPTY_QQ_REVIEW` の actual RAW24 全域全射と
`RAW24_PRODUCT41_BRIDGE_REVIEW` の点別証明に hash で戻した。
その全域は P=ad、B0=aB と変換すれば上の FULLU と同じである。
橋渡しに使う六式の変換 det80a、ℓ 座標の係数 (-16CJ/243)、
四つの modZW pivot (5a,aC/240,16C^2/81,160aC/81) はいずれも R-unit。
各点の restriction kernel の17方向と RAW24 の24方向から下界41、
係数10aの非零 syzygy から上界41となる。
これは abstract free row module の階数や generic rank だけの主張ではない。

さらに実際の `code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json`
の最後の関係を読み直した。非零行番号は

`[7,8,10,12,14,15,19,21,22,24,26,31]`

のみで、全て元の42 product 内にある。削除行14の係数は正確に12a/5。
従って選択41 product も元の全42 product と同じ像を張る。
genuine 行を利用した product の置換ではない。
元63行の六関係の unit triangular minor (-12a^3/5) も既受理の全91列恒等式に固定し、
上界57を保った。

## φ の product 像と Cdelta

元 physical91 の順序と full57→63 の全正規化は、独立受理済み
`FULL57_SOURCE_LIFT_REVIEW` に固定した。
Z=W=0 の後の L 座標の逆分母は (-40C) だけであり、実際の式は

\[
q_1=L^2/8100,\qquad q_3=L(eL-400ax_0)/8100.
\]

これにより ((L^2,Lx_0)) の degree3 の九単項式を除く。
80a の元六式恒等式から q2、q4 も ZWを法として(q1,q3)に入り、
全 product の像は q0 と scaled q5 の五 light 倍 Q10 の span に入る。
ここで19minorなどは反転しない。残る pure26 と H5、および

\[
\psi=-w_\delta HZ+z_\delta HW
\]

を合わせた写像を φ とする。products は H を含まず、ψ 像は0である。

全16 genuine の clearing (-64000C^{3-\mathrm{cancel}_i}) は R-unit。
relative1・6では同じ δ=-64000C²であり、実際の

\[
C_\delta=\delta(\text{full57 row47}-\text{full57 row42})
\]

は φ の核にある。元の行が零という意味ではない。
literal source の (c_Zz_\delta+c_Ww_\delta=\eta)、
η=192aFB²C³/125 は R-unit である。

raw7 cleanup の符号も、新しい自由可換多項式による小検査で確認した。
(b=c_ZHZ+c_WHW) とすると

\[
\eta HZ-z_\delta b=-c_W\psi,\qquad
\eta HW-w_\delta b=c_Z\psi.
\]

従って H5 と ψ が0なら ((b/\eta)C_\delta) を引いて全raw7 Hを消せる。
純三次 tail も同じ実際の source 行の差として運ぶので、省略していない。

任意の FULLU 幾何点で old25 の25行が独立なら、保持 genuine15 の product を法とする
関係は φ で Q10との関係になり、全 genuine 係数が0となる。
Cdeltaを加えた関係にも同じ議論を適用すると、その係数だけが残り、η-unit の H pairから
これも0となる。従って genuine16 は product を法として独立であり、
実際の product41 と上界57を合わせて actual full57 が従う。

## 8＋2個の global unit 操作

保存された八 pivot の分子・分母を、それぞれの named-unit 因子と独立な小多項式積で
照合した。pivot は順に

\[
-5/4,-5/4,5a,-5/4,25a/4,-a/6,-a/6,2FB^2C^2/15625.
\]

各 forward は pivot行の unit除算後に保存 pre-coefficient を使う row shear、
inverse は shear の符号を反転してから元の unitを掛ける操作である。
この一般形の両側逆も Laurent多項式の2×2模型で確認した。
checkpoint の選択8×8 block の対角・下三角36セルだけを別 parser で読み、
対角1・下三角0を確認した。上三角逆は strict partの有限和であり、新しい分母はない。
全 checkpoint 係数と保存八step自体は、既存の独立全係数検算に hash で固定した。

GL6、GL3 の実際の forward/inverse は、それぞれ既受理の独立 exact peer を参照した上で、
source の (I+u e_0^t)、(u_0=0) と inverse (I-u e_0^t) の形を直接確認した。
target determinantは -1、-9 の非零有理数であり、全新規 leaf 分母は QQ と a の高々二乗だけ。
3914 DAGの全命令は和・差・積または凍結 source leafで、逆演算命令を持たない。
2138-node prefix、old/new source・target順、15×22の全330 root を確認した。
両方向で heavy列 `[27,28,29,30,31]` は literal identityのままである。

Q9の選択九pure係数の外にある八tail `[1,2,3,4,6,7,8,9]` も全て残る。
tail2はその一つであり、八tailを零と仮定していない。
H3の `[28,30,29]` と full57/63 source cleanupも維持される。
従って **FULLUの全剰余体で**

\[
\operatorname{rank}(M_{25})=10+\operatorname{rank}(M_{15})
\]

が成立する。

## 独立係数証明との結合

DIRECT peer `f946ca63…` の環Rでの恒等式は

\[
A D_{15}^{norm}=-T D_{12}
\]

であり、Aを先に逆にしていない。
ここで T は `KAPPA_EXTRA2/root239.json` の実際の7208項多項式であり、
Q2補助式900a−126e−121ではない。独立受理済みの T 非零性と、今回の一係数 peer
`c5ebf88f…` による実際の D12 非零性から、整域Rの **非空な開 (D(TD_{12}))** を得る。
そこでは右辺が unitなので A と (D_{15}^{norm}) も unitとなり、current15 rank15、
old25 rank25、actual full57 が順に従う。追加の κ3・A・P・旧Ω因子の逆数は不要。

旧Ω条件付きファイルは変更していない。新開が幾何学的に「真に大きい」こと、
actual bad Fitting idealとの等号、全FULLU57、全quartic239、exhaustionは主張しない。

実装時の最初の試行は `remaining26_order` の説明文二個を26個の index と誤読した assertionで
停止した。試行コードと failureを `_ATTEMPT1` に保存し、正しい
`remaining26_full35_column_indices` へ直してから全検査が通った。
元の source・PRIMARY の bytes は変更していない。
