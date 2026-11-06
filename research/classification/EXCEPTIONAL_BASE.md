# 特殊一次部分の被覆準備と、一般五点開に必要な底因子

2026-11-06 / v64-classification-first。課題E0。これは全被覆Eの完了ではない。

## 対象と今回の確認範囲

元S61を $X=V(I)\subset W\times\mathbb A^{45}$、$W=\mathbb A^{16}$、
射影を $\pi$ とする。$I$ は元134式を全て含む。
重み1の歴史的u座標は

`u4,u5,u12,u13,u14,u25,u26,u27,u34,u35,u36,u37,u48,u49,u50,u51`

である。不変式G28で用いるprimitive betaの規約とは、既存の可逆な線形adapterを通して比較する。

[52件の台帳](exceptional_base_inventory.json)は、凍結された51件と追加NEW25の全圧縮hash・復号hash、
16出力、全逆元節点とその源変数依存、像次元の証拠を今回照合したもの。
像次元の証明そのものの再計算や全幾何の独立査読は行っていない。
元134式のbeta零制限だけは今回正確な有理係数で構成した。
新しいrank、factor、minor、jet計算は使わない。

## 既知像の区別

既知52成分のうち、$W$ を優越するのはT01、T39、T40の三件。
残り49件の閉像 $Z_i=\overline{\pi(C_i)}$ はproperである。
次元別件数は $16:3,15:7,14:28,13:4,12:1,11:4,10:3,9:2$。
これらは既存の関数体rankと幾何上界の証明に依存し、今回の有限標本からの推論ではない。

NEW25の閉像の既約方程式は既にG28として得られている。
古いdense-map recordの `firstjet_image_hypersurface_equation_constructed=false` はこの後続結果を反映していない。
既存の完了済み証明を優先する。
G28は78項の不変式表示で、次数28、複素数体上既約。
証拠は
[既存証明](../notes/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63.md)、
[独立査読](../notes/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63_REVIEW.md)、
`reports/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63/equation.json`。
これはactual42 scalarのdelta_E次数28とは異なる多項式である。

NEW25の像には $\Delta\ne0$ の既知点がある。
従って特殊成分の底を $V(\Delta)$ だけとみなすことはできない。
この事実は既存のdense-map独立査読で確認され、台帳に束縛している。

## 一般五点開の分母に対する必要条件

**補題。** 有限型被約スキーム $X$ の相異なる既約成分を $C_i$ とし、
$D=C_{T01}\cup C_{T39}\cup C_{T40}$ とする。
非零の底多項式 $f\in\mathbb C[W]$ が

\[
 \pi^{-1}(D(f))\subset D
\]

を満たすなら、$D$ に含まれない各成分 $C_i$ に対して
$f\in I(Z_i)$、$Z_i=\overline{\pi(C_i)}$ である。

証明。$f\notin I(Z_i)$ なら $C_i\cap\pi^{-1}(D(f))$ は非空開集合。
他方、相異なる既約成分の有限和に $C_i$ は含まれないから、
$C_i\setminus D$ も非空開集合である。
既約性によりこの二開集合の交わりは非空で、仮定の包含に矛盾する。□

**適用。** 既知52成分の成分性・相互識別とS61での対応を既存証明から用いると、
全S61の既知五点を保証する底だけの分母は必ず

\[
 f\in\bigcap_{i\in\mathrm{known49}}I(Z_i)
\]

を満たす。特に、primitive betaの規約で **$G_{28}\mid f$**。
これは新しい分母候補の整合性検査であり、$G_{28}$ を反転すれば未知枝が全て消えるという十分条件ではない。
一般関数体だけでの証明にはこの因子を逆元にできるが、分母を払って底の開へ移す時には
$V(G_{28})$ の被覆義務を残さなければならない。
一般五点の仮定を、新しい成分の存在で反証したとの主張でもない。

## 49個の閉像を有限な方程式問題にする

源パラメータをx、各演算節点に対応する変数をtとする。
凍結DAGの `q(c),v(j),+(J),*(J),inv(j)` に、それぞれ

\[
 t_i-c,\quad t_i-x_j,\quad t_i-\sum_{j\in J}t_j,
 \quad t_i-\prod_{j\in J}t_j,\quad t_it_j-1
\]

を対応させ、16出力に $\beta_k-t_{o_k}$ を加えたidealを $\Gamma_i$ とする。
すると **$I(Z_i)=\Gamma_i\cap\mathbb C[\beta]$** である。

証明。節点順に消すとgraph環は源の多項式環を指定逆元で局所化した環に同型である。
源の既存の非空dense開上での写像に一致するため、betaへの準同型の核が閉像のidealになる。
beta出力に不要な節点を除いても、定義域のdense開を大きくするだけなので閉像は変わらない。
これは源の既約性と各除外元が零関数でないことを使う。
DPMでは係数体 $\mathbb Q(\sqrt3)$ の複素埋込みを一つ固定する。
$r$ を超越変数としたり、$\mathbb C[r]/(r^2-3)$ をdomainとみなしたりしない。□

実装もこの特殊化を式へ反映する。DPMのarchiveが源入力として持つrは、
export時には源変数一覧から除き、対応節点を `coeff(rho)`、すなわち
係数体 $\mathbb Q(\rho)$、$\rho^2=3$ 上の式 $t_i-\rho=0$ に置き換える。
他の源変数の添字も付け直す。両共役を複素数体上の一つのdomainとして扱わない。

この有限idealの表示は、像の展開方程式やeliminationの完了を意味しない。
NEW25の例では全6079節点のうちbeta出力の祖先は1350節点であり、
全係数とinverseを保ったgraph recipeを今回保存した。
全52件を同じ実装で取り出せる。
G28の既知表示があるNEW25で、これを再implicitizationする必要はない。

## 源の分母と、底・ファイバーの補集合

catalogueのvalue pivot、value Pfaffian、Sp4大セル行列式、period分母、scaleは源の関数である。
それぞれの零集合は源のchartの補集合で、Wの多項式が明示されたわけではない。
これをそのまま $V(f(\beta))$ に置き換えない。
源のdense-map theoremは既知成分の閉包を同定するが、源のchartから失われた特殊点や、
その閉像の上にある全Xのファイバーを自動的に覆わない。

またGの途中で未知qを逆元にした場合、その補集合はX内部の閉集合である。
底の補集合として扱うには、その射影像・閉包または別の全点被覆を証明する必要がある。
FULLUのa,d,F,B,C,Jも特定無限遠表示のパラメータであり、直接Wの除外因子とはしない。

## 実際に残す閉集合のキュー

1. Gの証明が与える底だけの $f_G$ の $V(f_G)$。現段階では $f_G$ 自体が未確定。
   上の補題から49像、特に $V(G_{28})$ を含む必要がある。
2. $V(G_{28})$ の全Xファイバー。既知NEW25があることと、その上を全て覆うことは別。
   残り48像のprime contractionも台帳から取り出せるが、未知特殊成分をこの49像の和に限定しない。
3. $\beta=0$。今回保存した元134行の制限は45変数の二次70行と零64行。
   旧 `Q45_22` はこの層の被約22成分の尽くしを収録するが、22個の全族から
   既知52成分への埋込み・閉包対応を主定理で使うには対応表を監査する。
4. $V(\Delta)$ と、Gから生じる他の底因子との交わり。既知T01の退化だけで全Xを尽くさない。
5. 未知qを反転した全chartの残部。源chartの分母との混同を避け、底への像がproperかを判定する。

各層で元134式を底idealでbase changeし、既知族閉包の和への包含を示す。
源へ引き戻した式の零だけでは、その底上にある未知点を排除できない。
層の一般部分を覆えたら不成立集合を真の閉部分へ送り、次元減少を記録する。
これはEの有限再帰の準備であり、Eの全停止を保証する証明ではない。

## 再現

新しい出力先を指定する。全て既存の凍結入力を読む。

```sh
.venv/bin/python code/classification_exceptional_base.py inventory --output work/E0-inventory.json
.venv/bin/python code/classification_exceptional_base.py beta-zero --output work/E0-beta-zero.json
.venv/bin/python code/classification_exceptional_base.py image-graph --component R4_POSITIVE_WEIGHT6_25 --output work/E0-graph.json
```

保存した制限は `reports/CLASSIFICATION_E0_BETA_ZERO_20261106.json`、
graphは `reports/CLASSIFICATION_E0_NEW25_IMAGE_GRAPH_20261106.json`。
今回の構造・係数検査と、既存証明を依存として使う部分を区別する。
