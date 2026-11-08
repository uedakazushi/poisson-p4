# Quadric cone・rank二・C零の全空性：独立純証明レビュー

2026-11-08。PRIMARY E。**C零枝だけをPASS**とする。
非零semisimple／nilpotent Cの全空性は本レビューで受理しない。
新しいCAS、operator、kernel、minor、point、profile又はjetは実行していない。

## 1. 全入口と正しい先行依存

`QUADRIC_CONE_RANK_TWO_ZERO_C_ACTUAL_FACTOR.md`（15bb7516…）を独立読解した。
全元33源・vertex・Euler商・actual normalizer・Hamiltonian pencil消去は、既別作者
`QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md`（4c1f15c8…）の全源証明を使う。
そこで元JacobiからDpencil零を得た**後**に限り、strict h一chartで

    P=B(a,b,c)+(Y(a,b,c)+vT(a,b,c))∧∂v

となり、degB≤2、degY≤3、degT≤1である。
入口でBのv非依存性を仮定した旧失敗sourceは使わない。

C零／Fv零では、先行独立全126 actual係数証明からh因子がある。
Fv非零はSO(ac−b²)とactual scalar／normal scalingの全逆により、
q-dual norm非零のb型、nonzero isotropicのc型へ尽くされる。
零Fvは除外して消すのでなく、上のh因子枝へ戻す。
restriction零divisor Rのreduced性又はvertex一般位置を仮定しない。

## 2. 原33方向を読むと、C零では本当に全six base sourceになる

`CLASSIFICATION_QUADRIC_CONE_RANK_TWO_PENCIL_20261107.json` のcases[2]、
`pencil_B0_point_and_all_directions` 全34個（index零はparticular）を確認した。
全Bに寄与する自由係数は正確に

    A=p1, d=p2, D=p8, r=p7, s=p14, t=p20

及びFvのp28,p31,p33だけである。
**C零caseのFu-only directions p26,p27,p29,p30,p32は全Bで零**であり、
非零Cで現れる追加base方向を無断で消しているわけではない。
他の全方向はY,Tへ残す。

Fv=bではp31一、p28=p33零で、全Bは

    Bab=A−rb+sa, Bac=d−rc+ta−2q, Bbc=D−sc+tb。

Fv=cではp33一、p28=p31零で、全Bは

    Bab=A−rb+sa+2q, Bac=d−rc+ta, Bbc=D−sc+tb。

q=ac−b²。ここでd、D等はfree constantであり、正常化したconstant bivector Cとは別物である。
全Y,T、全sigma及びその共鳴を一つも零へ制限しない。

## 3. 元base Jacobiの全係数と、全actual factor

bvec=(Bbc,−Bac,Bab)とする。
Fv=bの`bvec·curl(bvec)/2`を手で展開すると

    −tab−rbc+2sb²−Da+Ac+sd−rD−tA。

ab,bc,b²からt=r=s零、a,cからD=A零。
dにはunit条件を付けず、全Poisson sourceでB=(d−2q)∂a∧∂cとなる。
Fv=cでは`−bvec·curl(bvec)/2`が

    rc²+tac−2sbc+2Db−dc+rD−sd+tA。

c²,ac,bcからr=t=s零、b,cからD=d零。
Aを零又は非零へ限定せず、B=(A+2q)∂a∧∂bとなる。
この二つの式の外積符号を独立に確認した。

全四変数Pfaffianは前者で符号を除き
`(d−2q)(Yb+vTb)`、後者で`(A+2q)(Yc+vTc)`である。
全Y,Tのままの恒等式なので、higher compatibilityの追加条件で積が壊れることはない。
第一因子は全d,Aで正確にdegree二、cofactorは高々degree三。
元degree五へのhomogenizationは`d h²−2q`／`A h²+2q`というproper quadric因子を与える。
cofactor零ならactual D零でrank四の対象外である。
Euler gaugeはactual E∧P²を変えず、全actual normalizer逆も非定数因子を保つ。
従ってrank二quadric cone、C零の**全rank四・被約既約quintic源は空**。

## 4. 非零Cの別の入口欠陥は混ぜない

元PENCILのSS caseには、例えばp26のBbc=4aとp32のBab=−4cがある。
nil caseにもp26,p27,p29,p30の追加linear base方向がある。
従ってそれらを消す全域Jacobi帰結又は実正常化の証明無しに、非零Cを上のsix型へ移せない。
authorもこの欠落を確認し、旧six-base及び下流全closureを不採用として別訂正を準備している。
このことはC零caseで上記方向が実際に零である事実に影響しない。
今回のreceiptはC零だけを狭く認証し、非零Cの全source、既約quintic不存在、
全quadric-cone E完了又は新成分数を主張しない。
