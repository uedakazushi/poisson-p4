# m4 repeated leading・全 generic211 source の G 上界

2026-11-07、主義務 G。producerによる全幾何還元・exact Fraction算術。
実 D は既約被約 quintic、S は restriction rank 二の plane、zero cubic R は reduced reducible とする。
S沿い位数四、leading binary quartic の generic root partition が 2+1+1 と仮定する。
**この全 source の projective PGL 像は次元高々23、特に G 非優越である。**
image-rank 一は既essential-three八成分へ、rank二・三は下記六型へ還元する。
generic22/4も必ずimage-rank一なので同じ既成分被覆を持つ。
generic31 image-rank二、m2/m3、非被約R、他支持と全Eは明示残余である。
この結論の別作者独立確認は別receiptを要する。

## 1. 全幾何と接続の適用範囲

`PLANE_M4_REPEATED_LEADING_GEOMETRY.md` §§1–5 の elementary P1 Bertini、
fixed square factor、二つの collision divisors、regular symmetric2 pencil の distinct/Jordan 両型を使う。
normal quartic に三つの distinct roots があるので twisted GL2 Lie stabilizer は零。
全 parallel connection は同note§7の式で一意に復元され、nonlog 接続も含む。
rank三 universal net の collision conic と marked line は tangent であり、
transverse conic+line の高rank枝を誤って追加しない。

実行では fixed double root をuとする convention に全型を揃えた。
Jordan型のG=xu²v²+yu³vは幾何noteの表示からu,vを交換したもの、
三used-lines型のG=u²[xu²−2xuv+(x−y)v²]は旧normal座標を(u−v,v)としたもの。
いずれも全normal GL変換であり、sourceを減らさない。
最初の四型の全parallel-leading kernelは二次元、後二型は一次元で、全非零振幅を保持する。

## 2. 一回の確定全六型計算

有限契約 `jobs/plane_m4_repeated_211.json` の gate後、reusable runnerの
isolated `--repeated-211` を実行した。
最終report `reports/CLASSIFICATION_PLANE_M4_REPEATED_211_20261107.json` は
4.712815秒、RSS37,504KiB、全算術PASS。
前のconcurrent実行sourceの全七FunctionDef ASTは不変で、新しい三関数とdispatchだけを追記した。
旧producer/report/証明書、essential-three全八packetは一切上書きしていない。

最初の実装実行はJordan型でf/Tの多項式性を誤仮定して停止し、final reportは書かなかった。
f*dT/Tを微分後に約分するだけの修正を行い、新source hashとgateで全六型を再実行した。
失敗実行を証明とは呼ばず、以下は確定reportだけを使う。

各型について以下を全て保持した。

- fβの全polynomial one-form、全12 curvature numerator zero。
- H0Ω1(3)とtracefree linear fieldsの完全18×8写像と四entryの逆解。
- actual D2=0によるraw normal-normal quadratic三係数の一意解。
- 全150 P0、全350 Jacobi、actual Pfaffian zero。
- 全十五leading coefficientsの完全parallel核とrank minors。
- 全三十weight一source、全九weight二source、全374×30と全350+shape×9行列・核・minor。
- 全P1二次項とoriginal-row compatibility functionals。

raw normal-normal bracketを初めから零とはしていない。
triangle二型はそれぞれuv/4、−uv/2、tangent-conic netは3u²/4であり、
concurrent三型だけが零となる。

## 3. 三つの全次元上界

全weight一JacobiとD3=0の必要核をK1とする。
任意P1∈K1上でweight二式は同じ定数operator L による L(P2)+Q(P1)=0。
各fiberは空又はker Lのaffine translateなので、全parameterに一様な上界を得る。
追加high-weight Jacobiや実Dの既約性はsourceを狭めるだけである。

| 全source | dim K1 | rank L | P2 fiber上界 | fixed P0 source上界 | 実群fiber | 全projective像上界 |
|---|---:|---:|---:|---:|---:|---:|
| nonconcurrent diagonal | 6 | 8 | 1 | 7 | 10 | 21 |
| nonconcurrent Jordan | 6 | 7 | 2 | 8 | 10 | 22 |
| tangent conic net | 6 | 8 | 1 | 7 | 9 | 22 |

同じ固定linear operatorを使っており、generic rankを境界へ外挿しない。
この三型ではQ(P1)自体は非零だが、全compatibility残差は零であった。
必要上界に追加quadratic ideal、GBやjetは不要。

triangleの実群はbase torus二、normal diagonal二、normal-to-base shear六の十次元。
tangent conicではPGL2のmarked-point Borelが、symmetric-squareのbase actionとnormal actionへ同時に持ち上がる。
そのprojective imageはbase二＋normal scalar一の三次元で、shear六を加えて九次元。
base二とnormal actionを独立に数えない。
fullconnectionの一意性によりこの実群はfixedP0をprojectiveに保ち、全sourceを保つ。
PGL5×Xのfiberへの右作用はPGL5因子で自由なので、像は24+dimX−dim群以下。

## 4. concurrent三型は共通外積因子で全空

concurrent三型の粗source上界は12+6、粗像上界は31,31,32に留まる。
しかし、完全K1の全十二方向、P0及び完全kerLの全六方向は、
全てz-independentなfieldとの∂z外積である。
全original-row Q matrixも全係数零なので、任意P1でP2∈kerLとなる。
従って全必要sourceが

\[
 P=P_0+P_1+P_2=\partial_z\wedge W,\qquad \partial_zW=0
\]

を満たす。rank≤2、actual D=0によりrank四sourceは空。
この証明は各型の`entire_necessary_source_common_partial_z_factor=true`、
完全核・全quadratic matrix・P0全150係数から直接検査できる。
正確な粗像31/32を隠さず、不要なminor増加を行わない。

## 5. 結論と固定入力

rank一全源はessential-three theoremにより八projective23成分で被覆される。
rank二・三の六型は空又は像21/22以下なので、全generic211 sourceの像は≤23。
generic22/4もrank一へ既還元済み。
generic31高rank、m2/m3、nonreducedRは結論に含めない。

| artifact | SHA256 |
|---|---|
| code/classification_plane_reducible_cubic_m4_special.py | ca361245a1915735e31c5bbb8acd1170a0492b0ed73d6ca42f5fb58ae6c6b2fb |
| reports/CLASSIFICATION_PLANE_M4_REPEATED_211_20261107.json | 0e00f1eaeb31a0051dd641704c7b8c5e27f090de1ec3764bd0fd072dd7dc590a |
| research/classification/PLANE_M4_REPEATED_LEADING_GEOMETRY.md | e18c758206443d2f0c15916b4152466e52c480a69a50bf69b0ce7d99693a0f83 |

全source/接続/群の独立reviewは別receiptを保存する。新成分数は変更しない。
