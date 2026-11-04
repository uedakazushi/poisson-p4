# R434 を反転しない directP5/H5：独立監査

受理済み directP4/H5 の source と768操作の指定一成分 source cofactorだけを使う。新しい actual係数・minor・rankは計算しない。全74 typed nodeの型、外部hash、源/標的順序、block recipeと積順序を別実装で照合する。

## raw row replacement

保存3cofactorは実 current15行[5,2,4]に作用する。元のrhoを保つと、cuは文字どおり `(-81/800000)*delta2`。ほか二係数は pair adj2 の二つの mu に `−9/40000`、`−63/160000`を掛けたものとなる。raw compact pair行とu4行のmiddle三成分はliteral零で、pair二列は2対2 adjugate恒等式で消える。従って新row5は全heavy五列で零、pure17列は実Bhatそのもの。

raw row5の置換はhybrid rho scalingより前に行う。cuはD(Delta5)上単元。source行2,4は変わらないので、逆は `(newrow5−c1*newrow2−c3*newrow4)/cu`。独立な可換Laurent環で15対15の両逆積を確認する。

## P5 determinant と block inverse

P5はpure四行と新Bhat行、列[22,11,13,25,15]。detP5は普遍 bordered determinant式そのもので、Sauxとなる。この式は一度P4のgeneric可逆開で導けても、多項式恒等式として整域全体へ降りるため、D4/R434の逆元を導入しない。新しいAdj5も polynomial operatorのまま保持する。保存された7022項のQQ係数を本監査で再展開したとは言わない。actual source SLPの非零性には既受理の独立 symbolic-modp 証明を使う。

P5とH5の block K10は `[[P5,0],[C55,H5]]`。H5の既存逆の全依存closure、およびprefix/source25輸送のinverseを調べ、FULLU既知単元とDelta5だけを使い、R434/D4/h22の逆を含まないことを確認する。自由な非可換block代数で、P5/その逆、H5/その逆の関係だけからK10の左右両逆とSchurの全block等式を検証する。

元15対22→actual同cokernel5対12は厳密に `FULLU D(Delta5*Saux)` 上で成立する。R434=0を除外する必要はないが、その側との交叉の非空性は別の問題である。残るsourceはcurrent15[8,9,11,12,14]、original25[18,19,21,22,24]、targetはoriginal32[1,2,3,4,6,7,8,9,14,17,18,19]。八つの追加Q9pure tailは全て残る。

rank5・全FULLU57・quartic239は未計算。Delta5=0、Saux=0と他のheavy chartは保持する。R434のfactorについてこの監査は新しい主張をしない。
