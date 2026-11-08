# R=2L+M・m三SF：最後のc-dependent源にも実z因子

2026-11-08、`/root/quartic_review`。追加の演算job、M1、点、P2解法、jet又はGBは実行しない。
既独立に検証した完全kernelの三つの零rowだけを使い、C=a²,g=2の**全既約quintic枝は空**と証明する。

入力は `CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_REPLAY_20261107.json`
の `cases[15]`, type=`C_DEPENDENT_A2_G2`。
full sourceと全kernelの独立算術は
`CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SEPARATE_REVIEW_20261107.json` 及び
`PLANE_DOUBLE_LINE_M3_SEPARATE_REVIEW.md`で既に受理されている。
本稿は当時のG≤24の上界からEを推測せず、保存済み全源に対する新しい因子含意を用いる。

## 全normal rowの三零係数

座標は(a,b,c,z,w)、normalはz,w。
全P0のz-rowは

    Paz=0, Pbz=0, Pcz=az, Pzw=2z².

元reportのfull_source1のパラメータ順は

    l0,...,l8,n0,...,n8,d0,d1,d2,b0,b1,b2,e0,e1,e2,c0,c1,c2.

この最後の三係数はそれぞれ

    P1_az=b0 wz+c0 w²,
    P1_bz=b1 wz+c1 w²,
    P1_cz=b2 wz+c2 w².

sourceにはP1_zwが無い。完全kernelは30×10のcolumn形式で保存され、
**rows27,28,29は全て零**である。
従って全必要K1上でc0=c1=c2=0であり、全P1のz-rowはzの倍である。
これはchosen vector又はleading imageの零ではなく、全10個のkernel generatorを同時に含む零rowである。
元350式とparallel条件に対するそのkernelの完全性は既独立receiptのnonzero minor＋全row代入で証明済み。

全P2は元の9係数をそのまま許す。
P2は(a,b),(a,c),(b,c)の三base pairだけに(z²,zw,w²)を持つので、
**affine particularを含めどのP2にもz-normal columnは無い**。
P2をhomogeneous ker J2=0だけへ狭める議論ではない。

従って任意の全Poisson solutionだけでなく、その必要な全線形supersource上でも

    P_zj ∈ (z)  (all j).

すなわちz=0は実tensorのPoisson hyperplaneである。

## actual Euler-Pfaffian

実五次式Dは五つのPfaffianをEuler係数で合わせたものとする。
z=0ではEulerのz係数を持つ項は零。他の四項のPfaffianは全てz indexを持ち、
そのrowが零になるのでやはり零。よってD∈(z)。
D≠0ならproper linear factorを持ち、D=0ならrank4の退化五次式にならない。
従ってこのwhole c-dependent源にはactual Dが既約quinticである点は存在しない。

squarefreeな**normal leading** cubicがあり得ることとは矛盾しない。
実際z times二次binary式はsquarefreeでも、全actual Dはz factorを持てる。
既leading-discriminantの非零を全actual irreducibilityへ読み替えない。

## 範囲と保存

本稿は旧19源の再計算を行わない純証明である。
新receiptは全bad三rowのliteral値、source monomial順、旧独立complete-kernel入力hashを保存する。
`C_DEPENDENT_A2_G2`の既約quintic枝を空に強化するだけで、
repeated cubic、m二、R=3L又は全分類の完了を主張しない。
