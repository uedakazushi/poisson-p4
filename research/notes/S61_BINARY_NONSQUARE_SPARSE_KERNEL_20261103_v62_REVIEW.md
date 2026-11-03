# V4 sparse kernel・三行削除・P=ad係数座標の独立監査

2026-11-03。凍結V4 packet、既存V3の原入力、および完成した
`S61_BINARY_NONSQUARE_CUBIC_REDUNDANCIES_20261103_v62` の公開証明書を監査する。
結論は、**指定開上の三つの厳密な行削除と、固定次数以下の核計算を受理**する、である。
全点でのcubic rank57、product rank41、全四次被覆、ambient除去はこのpacketの結論ではない。

## 有限行列とrank sandwich

左関係の各係数をa,e,dの全単項式（全次数≤D）で展開すると、
n*binomial(D+3,3)未知数の有限QQ係数行列C_Dになる。
core座標を有限体の点へ特殊化しているのではなく、全多項式係数を展開した有理数を法pへ落としている。
全分母がp単元であれば、独立なr本のmodp equationはQQ rank≥rを証明する。
一方、異なるfree coordinate上でidentityとなるN−r本のQQ核ベクトルが、
元の全91多項式列を厳密に零にすれば、QQ nullity≥N−rとなる。
この上下からrank=r、出力が**その固定次数以下の全核**と分かる。

rational reconstructionに高さの推測や一素数の一致だけを用いていない。
再構成候補が誤っていれば、全QQ恒等式の検査で拒否される。
early版も、未処理equationを含む元全列のQQ恒等式が上界を与えるため、全行走査前の停止が正当である。
full column rankが得られた負の結果は、有限体rankの下界だけでQQ核零を証明できる。
Mersenne素数2^61−1の積は128bitで計算し、一度のfoldと減算で正確に還元される。

## exact forced-zeroと独立試験

旧版のsingleton前処理は、QQ多項式を整列・同類項結合した後の**厳密な非零support**で行う。
ある係数equationに、既に零と分かった未知数以外が一つしか残らなければ、その係数はQQで非零なので未知数は零である。
この推論にmodpで消えた項を使っていない。
独立reviewerは次数4の90段、次数5の118段を元の全多項式から再構成して全段確認する。

新たに、5次元核を持つ小fixture、8係数全てがsingletonで消えるfixture、
法pでは零でもQQでは非零となる不運な素数の拒否を実行する。
実入力の次数2核は軽量に再計算し、既存の1次元核と一致することを確認する。
元の次数6・7やP座標の大きいrank計算は重複実行せず、凍結source・選択equation・log・exit0を監査する。

## 三つのQQ恒等式と開集合

原full63×91 matrixのhashは
`dae59904925260605628d925ae0565de084c9af0ae1bccf1774f1af067830ad6`。
各rowの正規化にはa,d,F,B,C,Jと非零QQ定数だけを使う。
ここでF=500a+90e²+3e−24、B=(50a+18e²−9e−2)d+60e−40、
C=100a+18e²−9e−2、J=5ad−6e−4である。
reviewerは保存した単元倍率について、原row=倍率×正規化rowを全成分でQQ再検査する。

原row48、14、55（0始まり）の係数は順に1、12a/5、1である。
関係の係数次数は順に2、3、6。後の関係は既に削除したrowを含まない。
三関係×全91列を、native discoveryと独立なPython Fraction多項式演算で零と確認する。
従って開a*d*F*B*C*J≠0の全点で三行を削除でき、残る60×91の局所化row moduleは元と等しい。
一般階数57への依存関係はなお三つ残り、60行の独立性は主張しない。

## 負の次数結果とP座標

元座標の**二行削除済みの固定61×91正規化matrix**では、次数≤4,5の追加核は零。
三行削除済みの固定60×91では次数≤7の追加核は零。
これらは全syzygy moduleが零という主張ではない。

P=a*old_dはa≠0上の可逆な係数座標である。
a^i e^j old_d^kはa^(i−k)e^j P^kへ移し、各rowへaの冪を掛けて分母を払う。
全entryの変換と記録hashを独立に再構成する。
新しい単元はa,F,a*old_B,C,5P−6e−4,Pであり、未確認の別因子を逆にしていない。
nativeの第三slot名dはこの箇所ではPを表す。
P座標での既知単元正規化も全rowの比例恒等式と因子を再検査する。

この**別の固定正規化matrix**に対する次数≤3,4,5,6の核は零であり、
有限未知数は順に1200,2100,3360,5040。
単元によるrow倍率やP座標への変更は、局所化row moduleや各点rankを保つが、
係数の次数制限付き核空間をそのまま保つわけではない。
従って各負の次数結果はそれぞれの行列・係数座標に限定する。

V4 READMEの「公開verifierは未実行draft」は保存時点の歴史的状態として正しい。
その後の完成した公開reportではcompileと全QQ再検証が成功している。
両状態を混同せず、凍結V4の101memberと公開reportの全hashを確認する。

独立reviewerは `code/s61_nonsquare_sparse_kernel_review_20261103_v62.py`、
receiptは `reports/S61_BINARY_NONSQUARE_SPARSE_KERNEL_20261103_v62_REVIEW/`。
