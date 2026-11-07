# 非被約R・m四高rank211の全Poisson源は空

2026-11-07、主義務G。producer全actual-D正常形証明とfull150 linear certificate。

**命題。** rank二plane S⊂SingD、D被約既約quintic、S沿い位数四、generic normal leading partition2+1+1、
leading coefficient image rank二以上、πSの零cubic Rが非被約、という全rank四Poisson sourceは空。
R=3Lの高rank枝はcollision支持だけで空、R=2L+Mの全枝は以下のcanonical sourceで空となる。

## 1. 全actual Dを二固定及び一wholelineへ戻す

`PLANE_M4_TWO_LINE_COLLISION_QUINTICS.md`の全幾何を使う。
rank三leading quadraticにはirreducible collision conicが必要であり非被約cubicのsupportに入らない。
rank二symmetric pencilの全semisimple/Jordan型及びmarked double rootを分類すると、
semisimple unaligned型は三used linesを必要とするため除かれる。
残るaligned semisimple及びJordan型は二相異used linesを必要とし、R=3Lにも入らない。
R=2L+Mの二つのordering、unused components、全nonzero amplitudesを保持する。

全actualDはPGL5と全体尺度を除いて、essential-four base(A,B,u,v)、vertex wを用い、

`sem Fε=u⁴A+u²v²B+εuv⁴+v⁵, ε=0又は1`、

`Jordan Fτ=u²v²A+u³vB+u⁵+τuv⁴+v⁵, τ∈C`

へ戻る。全base/normal shearsと逆を同noteに記した。
Jordanではirreducibilityに必要なu5/v5両端係数を一へ戻した後のτを、0/1へ誤って有限化しない。
全complex line、τ零及び全special fibersを以下で保持する。

## 2. canonical full150 sourceと唯一のparameter wall

任意のprojective quadratic Poisson tensorは一意なdivergencezero strict liftを持つ。
canonical Pfaffian–Casimir identityにより、actual D class [F]を持つその全源は

`LF={Π: 全五変数quadratic coefficients150, divΠ=0, Π#dF=0}`

に含まれる。任意Euler代表のdivergenceをglobal stabilizerとするshortcutは使わない。
この必要linear sourceにfull nonlinear JacobiとactualPf=Fを加えることは、源を狭めるだけである。

有限job`jobs/plane_m4_two_line_casimir.json`は、degree六Casimirの全1050係数とdivの全25係数、
すなわち全1075×150 matrixを直接作った。
`CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_20261107.json`には全basis、matrix、完全kernelと非零minorを保存する。
QQ(tau) generic minorの全factorはtau6だけであり、唯一のwall tau零も元の全matrixで再構成した。

| 全actualD | canonical LF dimension | 全whole-source次元上界だけでも十分な値 |
|---|---:|---:|
| sem ε零 | 2 | 19 |
| sem ε一 | 2 | 20 |
| Jordan genericτ | 3 | 21 |
| Jordan τ零 | 3 | 20 |

したがって単なるG≤24でも十分だが、完全kernelは次節のより弱いrank四空性を直接与える。

## 3. 全canonical sourceの共通vertex factor

full150 basisのうちbase四×base四blockの六pairに属する全90quadratic係数を調べる。
これはwを含むmonomialも全て含めた90であり、vertex dependenceを先に捨てたsourceではない。
上の各完全kernelの全columnで、それら全90成分がliteral零である。
generic Jordan QQ(tau) kernel及び唯一のwall tau零の完全kernelに対して同じことが成り立つ。

従って全τ及び両sem形で、canonical sourceの任意のΠは

`Π=∂w∧V(A,B,u,v,w)`

という共通constant vertex factorを持つ。
そのraw rankは高々二で、projective quotientのrankも高々二。
rank四projective Poisson tensor又は非零actual quinticを持つことは出来ない。
全kernelを使った結論なので、generic sourceを特殊fiberへ外挿していない。
actual Dの全normalization逆でrankとPfaffian零性は保たれるから、元の全sourceで命題が従う。

全90行のindexと各fullkernelのliteral零確認は
`CLASSIFICATION_PLANE_M4_TWO_LINE_VERTEX_FACTOR_CERTIFICATE_20261107.json`に別保存する。
これには新しいnonlinear elimination、minor又はjetの計算は無く、保存済み完全kernelからの直接の幾何的帰結だけがある。

## 4. 証明状態・reducedRとの接続・補集合

最終fullsource runは12.694秒、最大91,424KiBで全exact算術PASS。
実行producer bytesはreportへ埋め込み、以前の全二十八FunctionDef ASTを不変に保った。
matrix構成のzero-polynomial indexing及びJSON整数serializationだけを修正した二implementation failureも、
当時のsource/contract bytesとともにimmutable receiptへ保存し、失敗を数学の証明に数えていない。
独立reviewではcollision全支持、actualD全normalization、full150/1075、wholeτと唯一wall、
共通vertex factorをproducer import無しで再構成する。

本証明はactual D正常形を得た後にはRの被約性を使わない。
従ってreducedRのaligned sem/Jordan三profileと同じactualDを持つ全sourceにも適用出来る。
rank三conic型及びsem unaligned三line型は、既存reducedRの全六connection packetに残した別profileである。
image-rank一、generic22/4は既essential-three八成分、generic31は別の全-R二quintic theoremへ戻す。

全m四をまとめる際には、**nonredRのsquarefree leading枝を別に確認する必要がある**。
本命題だけでm四全R、m二/三、全rank二plane又は全G完了を宣言しない。
