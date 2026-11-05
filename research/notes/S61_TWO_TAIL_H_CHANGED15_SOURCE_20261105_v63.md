# 同じM12を使い三追加列を変えるgenuine十五行minor

2026-11-05。source/preparationのみ。新しいdeterminantや商環値は計算しない。
κ四係数が固定H_axis商環で消えたため、同じκのtail収縮を増やす方法は採らない。
既計算tail1・tail2とH27・anchorを再使用し、元十五行の別の選択minorを作る。

## 全sourceと正規化

旧M12のsource順を保ち、その後ろに残るpair三行を置く。

```
current15 row order = [0,1,3,6,2,11,12,9,7,13,10,14,4,5,8],
old12 target order  = [22,11,13,25,15,19,14,18,28,30,29,17],
new15 target order  = [22,11,13,25,15,19,14,18,28,30,29,17,27,1,2].
```

行変換は受理済みDIRECT_M12 peerの15×15 QQ[a,a^-1]行演算を使い、旧最後のcurrent14を
pair三行より前へ動かすだけである。具体的な最初の十二行は

```
Q9, (G18+2G39)/a, G39/a, G68, G20,
(-513/1600000)G115, (-2079/16000000)(G120-G99),
(-13851/88000000)G85, (-81/200000)G(current7),
(-21627*a/176000000)G(current13), (-2079/16000000)G99, G134.
```

残りはcurrent4,5,8そのもの。行名mnemonicではなくcurrent15添字と既受理operatorを正本とする。
受理済み行GLの逆も同じpermutationで運ぶ。旧順から三交換なのでrow determinantは符号反転。
source行GLはFULLUで可逆で、部分的なprivate cokernelを元の行空間と取り違えない。

H3 columns28,30,29と追加H27は1/L倍、追加tail1,tail2はB倍、他の列は1倍。
L=aFC²。従ってcolumn unit積はB²/L⁴。H31を使わずgamma inverseも不要。
全て完全反転積U6=adFBCJの既知単元である。
oldOmega、Delta407/Saux/N6/N7/κ/anchorや新しい不明因子を反転しない。

pair4の順はcurrent `[2,4,5,8]`。anchor列はliteral current15 column22である。
H27/L列は既存KAPPA sourceと同じく、compact14 rows `[1,3,4,7]` のraw H5 column0を
各rho constantで掛けた保存多項式。rhoのL因子をもう一回掛けない。
tail1/tail2列は既計算のB-cleared四列値で、別のa倍やκ倍は入れない。
この四列をP=[anchor22,H27/L,Btail1,Btail2]、A=P00とする。

## 符号付き多項式恒等式

pair四行は旧九pure列でrank一以下、H3三列はliteral零である。
最後のpair三行RiをA Ri-ci R_G20に替えると、旧M12の十二列は全て零になる。
行列式はA³倍される。最後の三行と三追加列の位置はいずれも `[12,13,14]`。
Laplace指数和は78で偶数なので符号は正。
残る3×3blockをAV-cuとすると、受理済み普遍恒等式

```
det(AV-cu)=A² det([[A,u],[c,V]])=A² det(P)
```

をそのまま使える。従ってsource domainでA²を消して

```
A*D15_new=+det(P)*D12.
```

となる。Aは保存156項の非零多項式だが、FULLUの単元とは仮定していない。
D15_newは元current15全行の選択minorに既知単元を掛けたもの。
GLOBAL全六因子source橋渡しによりD(det(P)*D12)はactual full cubic57の十分openである。
Q9の八tail、H3、Cdelta/rawH7持上げを一切落とさない。

## 固定有限商環での一回の証拠候補

既存軸e=d=0とprime32003、同じS=Fp[a]/monic(H_axis)だけを使う。
H27とanchorの保存係数、既計算tail1四列・tail2四列をaxisに制限した入力を用意する。
これは新しいphysical source評価ではない。aはSで単元で、全QQ分母も既存primeで単元。
Sを体・reducedと仮定しない。元の16成分・sourcehash・整数分母を保存する。

4×4det(P)を16subsetのsigned DPで求め、必ず結合積D12*det(P)を同じSで計算する。
別々の非零性だけではzero-divisorのため同じH成分上のopenが得られない。
結合積が非零で独立peerを得れば、H-squarefree/U6-coprimeからHのある二次元openで
このgenuine15minorとactual full57が成立する。H全体や全FULLUの被覆とは言わない。
結合余り零ならinconclusiveで止まり、H31への自動切換えや別tail/axis/primeを追加しない。
値実行にはfresh ROOT監査と一回GOを別途要求する。
