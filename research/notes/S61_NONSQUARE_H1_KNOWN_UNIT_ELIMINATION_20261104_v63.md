# 完全 H1 image の既知単元消去

前段の H² 単元消去後の actual free163→121 と同じ像を与える冗長190生成元
（quad57、cubic112、H²-free quartic21）を使う。H1射影で非零になり得るのは
q5H、全112cubic、全21quarticの134行で、残る56quad行はpureである。
純98係数は全てlazyのまま保持する。

凍結済みcompact恒等式から七つのsource行を
`delta*(G61-G20)*oldlight_j`へ変える。`delta=-64000*C²`はFULLUの単元であり、
この134行source変換のdeterminantはdelta^7である。そのH係数は
`(A_H*Z+B_H*W)*oldlight_j`で直接得る。117 target消去のH次数1 pivotは
L²、L*x0、ellZ、ellW、ell²である。この七行には最初、二番目、最後の係数が
なく、H射影はellZ、ellWの項を落とすだけである。compactの全source恒等式を
依存先として持ち、巨大な二行を別々に展開して相殺する必要はない。

各行に必要なscalar DAGだけをexact QQ(a,e,d)で評価し、最後に使い終わった
中間値を捨てる。ピボットはa,d,F,B,C,JとQQ定数の積と文字通り確認できる
係数に限定する。各完成ピボットごとに23係数、既存ピボットを引く係数、
元134行番号を保存する。全行に対するSchur操作はlazyで定義し、pure98部分も
元のfull sourceと同じ演算を受ける。未展開行・未完了ピボットを捨てない。

これは冗長source presentationの操作であり、n個のHピボットから自動的に
自由(163−n) presentationを主張しない。安全なlossless表示は
自由(190−n)→自由(121−n)である。既知kernel38、残り42は不変で、full121の
全点全射やpure98像の完了は主張しない。60秒、2GB、32MiBの一回の境界内で
得た完成段階だけをreportする。
