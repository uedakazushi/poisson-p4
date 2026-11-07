# normal degree六 rank零枝：別作者による部分還元の確認

2026-11-07。producer は root、reviewer は plane_collision_review。
`RANK_ZERO_NORMAL_SEXTIC_REDUCTION.md` を独立に読んだ。
**判定：記載された部分還元の射程で PASS。**
残る irregular ruled 解消の normal degree六モデルを不存在としたり、
normal degree六の全分類を得たとは解釈しない。
CAS、新たな係数 source、追加の数値上界は用いていない。

## 1. sectional genus四からの scheme-theoretic CI

normal S の特異点は有限個なので、一般 hyperplane 曲線 C はそれを避ける滑らかな
非退化 degree六 P3 曲線である。Castelnuovo 上界は g≤4。
g=4 では degree六の L=O_C(1) に対し Riemann–Roch と h0(L)≥4 から
L=K_C、h0(2L)=9 を得る。
既受理の integral-surface linear-normality により h0(O_S(1))=5。
hyperplane の完全列から h0(O_S(2))≤14 なので、非零 quadric が S を含む。
S の素性と非退化性は reducible quadric を除く。

rank4 quadric での class 制限を特に確認した。
Q は smooth quadric surface Q2 上の、vertex が一点の projective cone である。
vertex を除いた開集合は Q2 上の A1-bundle。
codimension三の一点を除いても Cl は変わらず、その開集合は smooth なので
Cl(Q)=Pic(Q\{vertex})=Pic(Q2)。
vertex を避ける一般 hyperplane Q2 はこの bundle の section であり、制限写像は
二つの ruling plane の class を Q2 の二 ruling line に送る。
従って C の bidegree が (a,b) なら a+b=6、(a−1)(b−1)=4。
これは a=b=3 を与え、S の class 自体が Cartier class 3H になる。

rank5 では Pic(Q)=ZH、rank3 では Cl(Q)=ZP, H=2P を使い、いずれも degree六から
[S]=3H を得る。後者二つの singular-quadric class group は
[Tironi, arXiv:1201.4508v1](https://arxiv.org/pdf/1201.4508v1) の
printed p.5 を直接開き、rank3/rank4 の記述を照合した。

normal Q 上の prime Weil divisor が Cartier class と線形同値なら、その divisor 自体が
Cartier である。従って S は O_Q(3) の section の零 scheme であり、cycle の一致だけではない。
O_P4(3)→O_Q(3) は global section 上で全射だから、S は scheme-theoretic CI(2,3)。
Koszul 完全列から q(S)=0 となる。

## 2. rank零 quintic の仮定が除く部分

既受理の全 Euler source 補題は、degree≥5、q=0、integral・非退化 S に適用できる。
その全源は I_S(2)⊗Λ²C5 である。
二つの独立 quadric があれば degree四の complete intersection に degree六の S が入り
矛盾するので、源は零又は Q*A のみとなる。
実五次式は零又は Q²linear となり、既約被約五次という仮定に反する。
従って今回の S は q(S)>0 であり、§1 により g≤3。
この矢印は smoothness、ACM、ordinary node、vertical divisor の仮定を使わない。

## 3. canonical Weil class と Schröer の正確な射程

C⊂Sreg なので curve adjunction は canonical **Weil** class に対して
K_S·H=2g−8<0 を与える。K_S の Q-Cartier 性は不要。
ample H は全非零 effective Weil curve に正なので、全 n>0 の nK_S は ineffective。
特に p_g(S)=h0(ω_S)=h2(O_S)=0 である。

[Schröer, *Normal del Pezzo surfaces containing a nonrational singularity*](https://www.math.uni-duesseldorf.de/~schroeer/publications_pdf/del_pezzo.pdf)
の revised 28 October 2000 PDF を独立に開き、§1、Theorem2.2 とその証明
（printed pp.2,4–5）を読んだ。
同定理は proper normal surface、perfect field、nonrational singularity、
全正 pluricanonical class の ineffective 性から q=0 を結論する。
del Pezzo 又は Picard rank一は同定理の仮定ではない。
著者の q は dim Pic0 だが、§1 は p_g=0 の下で Picard scheme が smooth となり
h1(O) に等しいことを明記する。今回 p_g=0 は上で確保済みである。
従って非有理特異点は q(S)>0 に矛盾する。

## 4. 解消に関する結論と残す枝

全特異点が rational なので、解消 f:S̃→S は q(S̃)=q(S)>0 を保つ。
f*H は nef and big で、projection formula により
K_S̃·f*H=K_S·H<0。
ここでは canonical class 自体の pullback を仮定していない。
全 pluricanonical class が ineffective なので κ(S̃)=−∞。
滑らかな複素曲面の分類から、S̃ は正種数 q(S) の曲線上の ruled surface の blowup である。
「rational」なのは S の**特異点**であり、S 自体が rational という結論ではない。

この irregular ruled 解消を持つ normal degree六モデルがなお未処理として残る。
滑らかな曲面の double-point formula を特異 S に適用していない。
Halanay の degree六論文本文は未読であり、その分類は今回のどの矢印にも使わない。
publisher PDF への既記録の403も、本レビューでは再試行又は迂回していない。
非normal S、他の degree、rank二 restriction、G 全体の完了も結論しない。
