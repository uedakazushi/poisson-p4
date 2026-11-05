# 固定順Bareissの全temporaryは1984bit以下

2026-11-05。旧一様2202bit上界を変更せず、固定row順から強い上界を別途証明する。
新しいmatrix値・axis値・determinantは計算しない。

L_iを各行の全entry全係数の絶対値和のsource上界とする。
任意minorの多項式全体のL1は、使う行のL_iの積以下である。
固定pivot段kで、既選択k行の積をP_k、pivot行をr、残りの行をiとすると、
Bareiss numeratorの両積のL1はどちらも

    W_(k,i)=P_k² L_r L_i

以下。従って引算前後の係数・多項式L1は2W以下。
quotient minorのL1はP_k L_r L_i、previous denominator minorのL1はP_k以下。
長除算の途中で差し引くpolyは、真のquotientの係数の部分集合とdenominatorの
部分係数の積なので、全累積subtractionのL1はW以下。途中remainderは3W以下である。
これはdivision回数を掛けた上界より強く、convolutionや部分longdivisionにも適用する。

sourceに保存された固定十二行順は0..11。各段で最大の3Wのbit数は
170,293,444,595,748,923,1130,1354,1566,1768,1984。
全division quotient自体は十二角minor上界1096bit以下で、入力scaleと出力もこれ以下。
よって全明示temporaryは1984bit以下、実装のhard guardは2048bitを維持できる。
新たな大きいbit capもSTOP-only completion仮定も不要である。

binary全witnessと全guard-state予算22.24MB、owned33275は旧binarygateのまま。
source atomのsigned軸評価、十二角determinantは依然として未実行。
