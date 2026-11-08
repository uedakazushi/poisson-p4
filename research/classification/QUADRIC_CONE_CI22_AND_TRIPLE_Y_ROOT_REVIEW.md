# 三つの全源因子定理：統合側の独立確認

2026-11-08。PRIMARY `E_reduced_nonnormal`。
対象を限った全源証明を受理する。全Eや総成分上界の証明ではない。

## 1. Quadric cone rank零：純手導出

既受理の全Koszul/Euler源・二つのA軌道・全Ore逆を前提とする。
`QUADRIC_CONE_ZERO_ACTUAL_FACTOR.md` のLie符号を独立に展開した。
半単純型で `Bxy=-xy, Byz=-yz, Bxz=d-xz-y²`、
`T=(-(c+1)x,-cy,(1-c)z)` とすると、y=0の全bad成分 `g=Yy|y=0` はaffine linearである。
固定cubicとradial quadraticのy成分はyで割れるから、全十五Y係数を保持してこの次数制約が成立する。
Lieのxyとyz成分は正確に

```
(d-xz)*g_z=-(c+2)*x*g,
(d-xz)*g_x=(c-2)*z*g.
```

`g=Ax+Bz+C` の係数比較で `(c+2)A=(c-1)A=0`、
`(c+1)B=(c-2)B=0`、`(c+2)C=(c-2)C=0`。
従って全c,dでg=0。cやdの反転、共鳴の削除は無い。

nil型では `Bxy=y²-px+d, Bxz=2yz, Byz=z(z+p)`、
`T=(k-cx+2y,p-cy+z,-cz)`。
z=0の全bad成分 `g=Yz|z=0=Ax+By+C` に対しxz成分は

```
(y²-px+d)*g_y=(k-cx+4y)*g.
```

y²、xy、yの係数から `B=4B, 0=4A-cB, 0=kB+4C`、従って全parameterでg=0。
元の追加互換条件を自由と仮定せず、必要一式だけで十分である。
両型でfull affine tensorの一row全体がそれぞれy,zで割れ、hyperplane上のrankは高々二。
全実Pfにその一次因子がある。chart w=1からhomogeneous Hへ戻すdense-open論理と、
元のλ零/C零のH因子補集合を確認した。したがって全quadric-cone rank零の既約quintic枝は空。
これは新算術ではなく、元Lie式の独立手導出と全源仮定の読解である。

## 2. Regular22_1 rank二CI22：全逆と疎Pf

`NONNORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md` の十源からの全逆を確認した。
二つの非零mapと二つの五源D零補集合を保持し、scheme上の二乗を一次式へ同一視しない。
`CI22_REGULAR22_RANKTWO_WHOLE_FACTOR.md` §2の六括弧を用いると、

```
F=P24*(x1*P03-x3*P01)+x4*(P01*P23+P03*P12).
```

自由rは最初から現れない。m²項は
`m²*x3²*((x1²+x3²)-(x3²-x1²)-2x1²)=0`、
sm項も `2x1²-q+(-x1²+x3²+x4²)=0` により消える。
s²項は `s²*q*(2x1²-x1²+x3²+x4²)=s²*q²`。
よって全自由r,s,mで `F=s²*x4*q²`。rank四点は既全H+2rank3Qへ戻る。
case Iと二つの零補集合の全D零は、保存配列の有限exact代入と既独立全逆を入力とする。
その小さい計算を全source分類そのものの独立再生成とは呼ばない。

## 3. Corrected-tripleの二Y壁：別Fraction算術

全九P1核・九P2 sourceの完全性は既独立receiptに依存する。
固定非零QQ九minorからP2の一意inhomogeneous解を戻す全消去と、反転無しの証明を読解した。
producerのSymPyと別に、既済標準ライブラリFraction engineから全十括弧を作り、
全350 original Jacobi係数及び全signed Euler Pfaffianを再構成した。
両壁とも全350が保存値と一致し、元row54/279は `2*p8²` 又は `2*p8²/3`。
従って全複素被約支持でp8=0となる。
実Dの全係数も一致し、そのp8零での全積式は

```
-z³*(-4a+p5*w+p7*z)*(-b+p3*z+p6*w)*(p0-4p1-p7)/4,
 z³*(p0-p2)*(-5a+p1*z+p5*w)²/5.
```

全振幅零を含め、D零又はproper因子を持つ。二つのwhole Y壁の既約quintic枝は空。
旧小切片、homogeneous P2零、chosenparameterの結果ではない。
事前gate後の独立算術は0.02355秒／RSS12,756KiB、全PASS。
`CLASSIFICATION_PLANE_TRIPLE_Y_WALLS_ROOT_REVIEW_20261108.json` にinline独立source、全入力SHAと検査範囲を凍結した。

## 4. Cusp-cone CI22：全rank四空性の追加読解

`CI22_CONE_CUSP_RANKTWO_WHOLE_EXCLUSION.md` の二つの全五sourceと元の零補集合を確認した。
case IIはrawで `V wedge partial_4`、case Iは表示された `h*P_h+sum(V_i wedge partial_4)`。
q1,q2への四Vの作用を直接微分して、二つは共通尺度 `2x2,-2x3`、二つは共に零と再導出した。
P_hのcontraction差がEuler倍という表示は、例えば第4成分で
`P_h#dq1=2x1*x3*x4, P_h#dq2=0` となることと整合する。
全五成分の保存多項式恒等式と元actualPf全零は有限certificateに依存する。
非定数q1/q2の微分をprojective tensorが消すので、全rank四開が空という論理を独立確認した。
rank四Poisson開で非零有理微分が消えないというchar0の事実を用い、
rawのEuler倍をstrict零と誤認しない。両openと全零補集合を戻すと全sourceがrank二以下。

## 5. 射程

四つの指定全源で新しい成分増分は零。
rank二quadric cone、他のCI22型、他のplane first-symbol、全Eの残源を消したとはしない。
個別因子定理を全source inventoryへ接続し、既約Dの行を空へ更新する。

## 6. Corrected-triple全点集約の読解

`PLANE_NONREDUCED_M3_GENERIC_TRIPLE_POINTWISE_COVERAGE.md` の全矢印を確認した。
正しいminus接続の200recordは12不整合、107全leading零又は実因子、36全special逆、43literal/newfixed包含、2Y壁に尽くされる。
全record indexの一意性と個数を独立に読み戻した。数表だけを数学の証明とは数えない。
36redirectはgeneric核の特殊化でなく全higher fiberの実normalscale同型であり、完全special核へ戻す。
43全包含は8constant-vector型、6strict-support型及び新fixedへ戻る。
constant-vector型は全actualDの二つの独立translation不変性からessential三となり、既全八成分定理を適用する。
6strict型の全18module/34壁は全補集合付き独立証明により因子・不整合又は新二成分に入る。
新fixedは元ideal内のD|z0恒等式、最後の二Y壁は本稿§3の別Fraction全源因子により空。
したがって被約既約quintic・非被約R・rank二plane・位数三・generictripleの全stratumは八essential三成分と新二23成分へ被覆される。
位数三でないleading零点は適切なm四以上へ戻し、そのrawsource全体の空性を主張しない。
旧plus90誤符号coverageは撤回のまま、全E一般へこの集約を外挿しない。
