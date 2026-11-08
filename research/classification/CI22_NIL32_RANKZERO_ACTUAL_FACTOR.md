# rank零 CI22・nil32 の全六源：二つの元Jacobi行で実二重quadricへ戻る

2026-11-08、PRIMARY E-FINITE-SURFACE。
**全sourceの必要因子定理を証明し、実既約被約quintic枝を空とする。**
同時に、この六源の全rank四Poisson点を既済の全 H+2rank3Q 支持へ戻す。
旧G次元上界だけを全点membershipへ読み替えず、新しい二行と全実Pfaffianを使う。

## 1. 受理済み全sourceとの正確な対応

`REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_REVIEW.md` §4–5 と
`CLASSIFICATION_REGULAR_PENCIL_CI22_REVIEW_V2_20261107.json` の Jordan partition `3+2` を用いる。
全Koszul/Euler源とabelian conormalの既存証明により、
rank零となるsingle-eigen nil32 CI22の全Poisson tensorは、実PGL変換後に

\[
 \Pi=q_1A+q_2B,\qquad
 A=a_0C_S+b_0C_R+c_0C_{R^3},\quad
 B=a_1C_S+b_1C_R+c_1C_{R^3}
\]

へ入る。六係数を一つも単位化しない。
S,R,R³とQの規約は凍結V2の全行列そのままである。
座標順序を x0,…,x4 とすると

\[
 q_1=2x_0x_4-2x_1x_3+x_2^2,\qquad q_2=2x_2x_4-x_3^2.
\]

ここで C_T=TQ^{-1}。非零上三角entryは正確に

\[
 C_S:(04)=1,(13)=-\tfrac12;\quad
 C_R:(03)=-1,(12)=1;\quad C_{R^3}:(01)=-1.
\]

従って f=a0 q1+a1 q2、g=b0 q1+b1 q2、k=c0 q1+c1 q2 とおけば、全六源は

\[
 P_{04}=f,\quad P_{13}=-f/2,\quad
 P_{03}=-g,\quad P_{12}=g,\quad P_{01}=-k               \tag{1}
\]

で尽くされ、他の上三角五entryは零。
この一つのcanonical pencilが全nil32型を代表するのは受理済みself-adjoint Jordan分類による。
全PGL飽和へ戻すときにも、以下の因子性とrank条件は保たれる。
旧Aut四次元のG計数は今回の全点証明には必要ない。

## 2. 二つの元projective Jacobi係数

single-cyclic convention の Jijk=Σ_l(Pil ∂lPjk+Pjl ∂lPki+Pkl ∂lPij) とする。
projective Poisson条件は E∧J=0。strict J=0とは置かない。
(1)から

\[
 J_{234}=0,\quad J_{134}=\tfrac12f\partial_0f,\quad
 J_{124}=-f\partial_0g,\quad
 J_{123}=-\tfrac12g\partial_1f+g\partial_0g+\tfrac12f\partial_1g.
\]

ゆえに四添字1234の全多項式は

\[
 (E\wedge J)_{1234}
 =-x_4\{(a_0x_2+b_0x_3)f+(a_0x_3+2b_0x_4)g\}.       \tag{2}
\]

q2はx0に依存しないので、(2)の次の**元の二係数**だけで十分である。

| 元350行、零始まり | 空間monomial | 全六parameterにおける係数 |
|---:|---|---|
| 310 | x0 x2 x4² | −2a0² |
| 314 | x0 x4³ | −4b0² |

従って複素Poisson閉点では必ず a0=b0=0。
generic parameterや選択点でなく全sourceに成り立つ。
全schemeでa0,b0が零、又はJacobi idealがreducedとは主張しない。
一次式のsyzygyやparameter GBを解く必要も無い。

## 3. 全実Pfaffianと零補集合

signed Euler Pfaffianを元の全十括弧から計算すると、全六係数について厳密に

\[
 F_P=-\frac{x_2}{2}f^2-x_3fg-x_4g^2.                 \tag{3}
\]

kは現れない。これはnormal symbolだけでなく、全actual quinticである。
a0=b0=0を戻すと

\[
 F_P=-q_2^2\ell,\qquad
 \ell=\frac{a_1^2}{2}x_2+a_1b_1x_3+b_1^2x_4.         \tag{4}
\]

q2はrank三で既約。ell非零は(a1,b1)≠(0,0)と同値。
従って全sourceのPoisson閉点は、次の二つの全補集合に分かれる。

* (a1,b1)≠(0,0)：F_P非零、一般rank四。実因子は正確に H_ell+2Q2。
* a1=b1=0：F_P=0、一般rank二以下。実際Π=k C_R³だけとなりdecomposable。

第一の場合、ellはQ2のvertex{x2=x3=x4=0}を含み、
dual relation `2*(a1²/2)*b1²−(a1*b1)²=0` より接hyperplaneである。
全a1,b1の零／非零を残し、a1又はb1を先に反転していない。

特に、当初の **actual Dが被約既約quinticであるnil32 rank零CI22枝は全空**。
remaining Jacobi互換性をさらに課す前の必要因子定理なので、その全互換性を
自由と仮定したり解けたと呼んだりする必要は無い。

## 4. Eへの全点の戻し方

rank四の全点は(4)のactual因子型 H+2rank3Q に入る。
これは既受理の全支持定理
`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §1、
`RANK_THREE_QUADRIC_SUPPORT.md`、primitive tangent／mixed／Tw0各全source theoremの対象である。
接配置・vertex包含・primitive／extra／Tw0の全補集合が既に被覆されているため、
この全六源のrank四点も名前付き既知成分の有限閉和へ入る。
normalizer subgroupの一orbitを仮定した新認識ではなく、**実全因子の定理への包含**である。
rank二以下の補集合はその既存の全支持義務へ戻し、rank四の新成分と数えない。

従って今回のE-finite-surface枝を閉じ、新componentの増分は零。
他のrank零scroll、nonnormal CI22のrank二source、quadric cone、
全E／全分類の完了までは主張しない。

## 5. 有限exact certificateと証拠の区別

`jobs/ci22_nil32_rankzero_actual_factor.json` のgate PASS後、
凍結V2のQ,N,S,R,R³だけから全150×6係数を再構成し、(1)と完全一致させた。
保存されたsource全体の350式を二十一parameter monomialへ制限し、
全350×21と全126×21 actualPf係数を保存した。
これは新しい150方向Jacobi operator、rank/minor探索、source profile、GB又はjetではない。
旧V2のq1CSについては、保存されていた元350係数と全て一致することも確認した。
(2)の全七十空間係数、二つの元行、(3)の全百二十六係数がexactに一致した。

算術は既存Fraction helperを変更せず利用し、**0.017718秒／RSS12,800KiB、全PASS**。
inline source、実行時helper bytes、全入力SHA、元row順序、全係数配列は
`CLASSIFICATION_CI22_NIL32_RANKZERO_ACTUAL_FACTOR_20261108.json`へ保存した。
本稿の手導出は統合担当も独立に再構成済みであり、正式別receiptは統合側で固定する。
gateやproducer自身の再実行を独立数学証明とは数えない。
旧review、旧runner、旧report及びcurrentjobを変更していない。
