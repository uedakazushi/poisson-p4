# R三重直線・m三squarefree leadingのrank二planeは空

2026-11-07、主義務G。producer全source証明とexact QQ leading-image certificate。

**命題。** actual退化因子Dが被約既約quintic、S⊂SingDがrank二plane、
πSの零cubicが三重直線3L、Sに沿うDの位数が三、generic normal leading binary cubicがsquarefree、
という全Poisson sourceは空である。
R非被約に沿う高位極を最初から許す。m二/四、m三の重根及びR=2L+Mは本命題の射程外。

## 1. 全first connectionを二つのactual正常形へ移す

`PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_FIRST_JETS.md`の全root-section証明を使う。
U=S\L=A2ではparallel binary cubicのdiscriminantがunit。
root incidenceはfinite étaleで三つのglobal sectionsへ分かれ、Pic(A2)零からpolynomial root vectorsを持つ。
pairwise determinantsがunitであり、base degreeの和が高々二なので、
三rootは全constant又はdegree pattern(0,1,1)の二場合に尽くされる。
後者はz、w+ell z、w+(ell+1)zとなり、ellはnonconstant affine-linear function。
actual base/normal linear changesでell=b/aへ直し、全mixed columnsをparallel tensorから一意に戻す。
twisted cubic actionのLie stabilizer零を使うため、nonlogな追加接続を捨てていない。

全first-symbol representativesは

`P0_scalar=3a² ∂b∧∂c`、

`P0_unipotent=3a² ∂b∧∂c+3az ∂c∧∂w`。

それぞれ全regular parallel leading spaceは四及び三次元であり、
leading cubicを一つのchosen pointへ固定していない。

## 2. 全actual leading mapはsquarefreeな項を持てない

global quadratic normal filtrationで全Π=P0+P1+P2を、P1三十、P2九の全sourceとして扱う。
全350 projective Jacobiのweight一と、D3が全parallel spaceに入るlinear条件を合わせた
complete matrixのkernel dimensionsはscalar二十四、unipotent十八である。
完全kernelをK1、元のactual Euler-Pfaffianの全24 D3係数へのlinear mapをL3とする。

`CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_LEADING_EXCLUSION_20261107.json`は、
full P0+P1+P2からactual Dを直接再構成し、全九P2変数がD3へ寄与しないことを確認した。
全350×30、全24×30 map、parallel-membership annihilator、完全kernelと非零minorを旧coarse reportと照合し、
全matrix product L3*K1を保存した。結果は正確に

| 全first connection | leading-image rank | 全actual D3 |
|---|---:|---|
| scalar | 0 | 0 |
| unipotent | 1 | λ a²z³、λは任意 |

である。unipotentの保存basisではλ=3(u15+u4+u9)。
全24係数のうちa²z³以外の行がliteral零であり、generic minorだけを用いた結果ではない。

scalarではm三自体に必要なD3非零に反する。
unipotentではD3零又は三重normal rootなので、generic squarefree cubicという仮定に反する。
任意のfull nonlinear Poisson solutionはこのlinear kernelに入るため、
weight二以上の全original Jacobiをさらに解く必要はない。
actual normalizationsの逆はsquarefree性とD3零性を保つので、元の全sourceで命題が従う。

## 3. 証明書・旧否定・独立確認の射程

reusable runnerのmodeは`--nonreduced-cubic-leading`、finite jobは
`jobs/plane_nonreduced_triple_line_m3_leading.json`。
60秒/1024MiBの契約をgate PASS後に実行し、1.721秒、最大79,640KiBで全exact検査PASS。
gateそのものを数学の証明とはしない。
全実行source bytesはreportへ埋め込み、以前の全二十七FunctionDef ASTを保持した。

旧`CLASSIFICATION_PLANE_NONREDUCED_TRIPLE_LINE_M3_SQUAREFREE_20261107.json`の
`all_G24_closed=false`とcoarse upper41/35は不変である。
大きなlinear sourceには必要なsquarefree actual leadingが存在しないという、より弱い十分条件で閉じた。
old reportを上書きしたり、同じdimension経路のminor/jetを増やした結果ではない。

別作者はfull350/30/9と全24 leading coefficientをproducer import無しのFractionで再構成し、
全root-cover/degree/actual normalizationの幾何と併せて監査する。
最終独立receiptの受理は、このproducer証明とは別に記録する。
全R=3L又は全非被約R/Gの完了は本命題から推論しない。
