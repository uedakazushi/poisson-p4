# 四monic quad の overlap5 unit と保存 qpair 入力の独立監査

2026-11-04。対象は `S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63`。
joint frame の実 Z/W補正を落とさず、固定98 pivot後の overlap18×141を監査する。
別に進行した qpair Gaussian は実行せず、本 packet が保存した **入力58×136**までを照合する。

## 独立 arithmetic と source normal form

全入力係数の denominatorは非零有理数×a,dの monomialである。
レビューは FracElement、author helper、author scriptを使わず、符号付き整数指数の
Laurent polynomial QQ[a^{±1},d^{±1},e]を Python Fractionの独立 dictionary arithmeticで計算する。
eの負指数は許さない。source numerator/denominatorから読み込む際にもこの条件を確認する。

actual六 prime quadと保存された monomial source definitionから全 product coefficientを作る。
normal formは leading-tail rewritingではなく、保存 unit square Uについて

    μU=row_C,
    μ_j=row_C,j−Σ_(i<j) μ_i U_ij,
    NF(row)=row_D−μRD

の sparse Schur式で構成する。98 pivotとdegree3の26 pivotは前の independent REVIEWで
全source/coordinate/leading coefficientを受理済みであり、今回も対角1と三角 supportを確認する。

全18 overlapのsource index、target141列、全2538 normal-form係数を再構成し、保存 inputと照合する。
同じdegree3 Schurで S=x0p1−Lp2のnormal formも再構成し、保存値がexactに非零であることを確認する。
単項式ideal M4による quotientや、Sを零にする仮定は使わない。

## 五つの既知 unit pivot とrank103

固定された五stepだけを読み、新しいpivotや単元候補を探さない。
各pivotは literal −80000*a/3で、FULLUの既知単元。
全行の正規化と保存された elimination coefficientsを照合し、18×141の after matrixと
18×18のsource row transform Tを独立に再構成する。
さらに、再構成後に **T*before=after の全2538係数**を別に計算する。
残る13 source行は全係数が零。五pivotが独立で、row operationは既知unitによる可逆操作なので、
overlap imageは全FULLUでrank5。その前のunit98 imageと合わせ、
**四prime quad の degree4 product imageは全FULLUでrank103**である。

これは五つの新しい image方向であり、新しい source-kernel方向を五つ得たという意味ではない。
十三の零行を元sourceへ持ち上げたpure-quadric syzygiesは、既存のwhole quadratic kernel22内にある。
既知38方向・残り42の数は変えない。

## q0/q5入力だけのsource binding

元q0,q5（prime p5,p6）の全29 cofactor monomialを同じ98 Schurへ送り、
上の五つの **既に固定された** image rowsでさらにreduceする。
全58×136係数を保存された remaining-pair-block.json と照合する。
その1286非零係数・62664 numerator termsを独立に数え直す。
保存されている21個のpositive known-unit entry witnessesも、指定座標と指定factor powersだけを検証する。
全候補を新しく探索したり、追加Gaussianを実行したりしない。

五unit pivotは、全 source221→target141表示にも同じtarget projectionとして作用させられる。
genuine112とquartic33は、そのexact projectionをlazyに作用させて保持する。
ただし、この REVIEW でその全係数、old38の新しいglobal complement、さらに小さいfree source basisを
構成したとは主張しない。特にfree183からimage pivot数を単純に引いて、新しい自由表示のrankを断定しない。

全四次141/136/239の全射性、全点cubic57、全kernel exhaustionは未証明。
四quadのrank103と保存qpair入力のliteral bindingが本監査の結論である。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_joint_Lell_overlap_unit_review_20261104_v63.py --output work/marked-square-v63/joint-overlap-review-REPLAY
```
