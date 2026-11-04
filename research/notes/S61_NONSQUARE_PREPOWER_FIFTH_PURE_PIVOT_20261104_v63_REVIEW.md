# H7 coverage と第五 pure pivot：独立監査

対象は凍結済み `S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63`。新しい actual係数・source query・rank/minor探索は行わない。全62 typed nodeを別実装で型・hash・index照合し、source/targetの変換を独立な Laurent 多項式環で検算する。

## H7 の外の幾何学的被覆

actual sourceから読んだ p,g1,g2,Theta と、原定義による B,C を QQ[a,e,d] で掛け直し、`750000p=g1*Theta` と `B²=Theta+C*g1*g2` を確認する。B は FULLUの単元なので g1 と Theta は互いに素。従って D(p)、V(g1)、V(Theta) が FULLU を覆い、p=0 の交叉片を落としていない。

D(p) の sparse pair determinant は `−6Bp²/5`。他の三列は七行 B にのみ支持され、pair二行では零なので、rank rawH は 4+rank B。g1=0 も既監査の別の単位 pair二行が middle三列で零であり、同じ分解が成立する。Theta 上の exact cokernel は既監査の `T/(P80P100)` で、B の選択 minor は単元である。従って rawH<7 の幾何学的支持は Bbad13 と Theta上二曲線で正確に尽くされる。

三集合上の full cubic57 は別々の独立 source/regular-representation/replacement 証明により全点で受理済み。よって未解決の実 full cubic 欠損は H7 開集合に限られる。Fitting idealのscheme等式は主張しない。

H7上で15 source→heavy5は全射だから、その kernelは局所自由10。Q9 sourceはheavy零かつ独立に受理済みの unimodular pure rowなので、sourceの直接因子と targetの単位方向を同時に除ける。残りは局所自由9→pure16となる。ただし projected full rankは実 cubic57 の十分条件に留める。

## 相対 6×13→5×12 の有限操作

受理済み direct P4/H5の E6×13 で v=E[0,9]を取る。vの source式は、実P4と二つの pair補正から作る bordered5 Saux に対し

`v=rho_u*Saux/(delta2*rhoLCM*D4)`。

D4はa⁴C³R434、delta2はDelta5の因子、rhoは既知単元。従って今回の chart は厳密に `FULLU D(R434*Delta5*Saux)`。h22は反転しない。Sauxはactual bordered-SLPとして解釈する。保存7022項の全QQ係数をこの監査で再展開したとは言わない。非零性は受理済みの独立 symbolic-modp source証明により標数零へ移る。

列9を先頭に移して E=[[v,r],[c,D]] と書く。独立な自由 Laurent 環で L=[[v^-1,0],[-cv^-1,I]]、U=[[1,-v^-1r],[0,I]] と逆を作り、全78成分の `LEU=diag(1,D-cv^-1r)`、左右両方の逆積を検証する。元source15/target22への埋込みと既存operatorとの積順序も別に照合する。source25への surviving5の係数式も保持する。

残り current15 sourceは[8,9,11,12,14]、original25は[18,19,21,22,24]。original32 targetは[1,2,3,4,6,7,8,9,14,17,18,19]であり、Q9の八つの追加pure tailは先頭八列に全て残る。

結論は指定chart上でのactual同cokernel5×12と `I6(E)=I5(T)`。Tの実係数やrank5、全域57/239は計算・主張しない。60個のbordered6支持上界は今回再計算せず、最小の既存見積180479が上限150000を超える停止状態を維持する。
