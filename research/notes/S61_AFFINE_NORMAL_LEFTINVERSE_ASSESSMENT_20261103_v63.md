# 線形normal行列のD0/D1 left inverseを再実行しない理由

2026-11-03。新しい129×42線形行列A(q)について、B(q)A(q)=I42が得られれば
元Iから全42normalを消せる。しかしD0/D1の同じ有限体計算を重複実行しない。

degree0のleft inverseなら、元134式のqに依存しない行結合で非零normal線形式を得る。
既存134二次最高部分のQQ rank134証明に反するので、**一般beta上でもD0は不可能**。

degree1のleft inverseならn=B(q)G(q)となり、Gは元134式の定数行結合だから、
元45変数のordinary degree3 Macaulay spanに42独立な線形式が入る。
しかし**同じ**beta中心・同じF32003で既存完全計算があり、
全6164行のrank6073、pivot counts=(0,0,135,5938)で、線形帰結は零だった。
第135式を加えたdegree3閉包もcomplete、rank6118、pivots=(0,0,135,5983)で線形帰結なし。
従ってその有限体fiberでD1は既に不可能と分かっている。

この有限体の負の結果だけを、一般QQ(beta)上のD1不存在へ延長していない。
係数線形系の全列rankをp単元minorで保証すれば一般化できる場合があるが、
その追加証明は今回は計算していない。
予測された負結果のための新たな120秒solverは**起動しなかった**。

receipt `reports/S61_AFFINE_NORMAL_LEFTINVERSE_ASSESSMENT_20261103_v63.json`。
既存報告と元入力/source hashを再確認した。古い大きなMacaulay計算を再実行していない。
全残余のunit I+(E42)=(1)、高次数cofactor、別の構造的な削減は依然別課題である。
