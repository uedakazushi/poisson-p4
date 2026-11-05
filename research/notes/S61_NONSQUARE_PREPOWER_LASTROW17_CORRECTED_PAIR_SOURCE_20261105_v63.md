# 列17の全 pair-tail 補正を含む四＋十 source 分解

2026-11-05。新しい係数は評価しない。非零 gate を満たした PIVOT14 の一対八を、
既存の最初の四角 block と残る十角 block に分け直す。

current15 の四行 [2,4,5,8] と original32 の四列 [27,31,22,2] を先に置く。
これは STAGED443 の A4 であり、known FULLU column normalization の下で
det A4=-T、A4逆=-adj(A4)/T が保存済みである。八個の fixed raw pure 列での
rank一を、残る列17に拡張して用いることはしない。列17の四成分全てを保持する。

残る pivot source 十行は [0,1,3,6,7,13,10,11,12,9]、最後の行は14。
残る pivot 十列は [11,13,25,15,28,30,29,19,14,18]、最後の対象列は17。
正規化した入力の該当十五対十五を

    [ A4  B   v4  ]
    [ C   D   v10 ]
    [ c   d   v1  ]

と書く。x=A4逆 v4、S10=D-C A4逆 B、b17=v10-C x、
d17=v1-c x、d10=d-c A4逆 B とする。求める一対八の列17は厳密に

    r17 = d17 - d10 S10逆 b17.

これで三つの追加 pair-row correction と G20 自身の効果を全て保持する。
M11 と G134 だけの十二角 border には置き換えない。

新たな十角逆を計算する必要はない。旧 PIVOT14 の L15,U22 と、既存の
column-normalization Q22 を用い、V=Q22逆 U22 と置く。
L の最初の十四行は元 source14 の成分を持たず、L[14,14]=1。
V の free 八行・pivot十四列は零、free八行・free八列は I8。
これらを有限な構造演算で全 slot 確認する。したがって selected normalized
十四角 K14 には

    K14逆 = V[selected targets,0..13] L[0..13,selected sources]

という既存 chart 上の両側逆がある。その右下十角 block が S10逆である。
PIVOT14 の元 source/target との block 三角恒等式から、r17 は元一対八の同じ
列17成分そのものになる。追加の anchor、g、p、R434、h22 の chart はない。

全 scalar node は既存 source/operator への参照または有限行列演算とする。
未知逆元を新しいアルゴリズム指示で導入しない。十角逆は上式で明示した既存
作用素の部分 block で、expanded adjugate や Gaussian を実行しない。
元25/full57/canonical63/pure84/raw7 は同じ last-source covector で持ち上がる。

列17の四 pair 成分の literal source は49 node/20 atom/2838 numerator項。
全十五成分は159 node/64 atom/10406項。これらは原 atom 数であり、x、S10逆、
r17 の展開量上界ではない。実 scalar job の残る gate は、既存 rational
coefficient を展開する schedule と全 lifetime/高さ/RAW の厳密な上界である。
本 packet は bounded typed construction だけを完了し、lastrow scalar worker の
実行を承認・実行したとはしない。全 FULLU57 と補集合は未決である。
