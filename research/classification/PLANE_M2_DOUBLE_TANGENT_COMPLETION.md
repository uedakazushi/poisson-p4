# rank二plane・tangent conic・m2 doubleの全G完了

2026-11-07。指定locusの全projective像は23次元以下。
個別成分数は増えず、全Poisson分類の完成を主張しない。

仮定は実D既約被約quintic、restriction rank二のPoisson plane、
reduced tangent conic+line zero cubic、normal位数二、generic double normal quadraticである。
primitive root度零又は一の全source、全parameter walls、全P1/P2を含む。

| 全leading模型 | 全sourceへの結論 | 証明 |
|---|---|---|
| constant h=f、off Z=0 | 像22以下 | 全30/9とactual filtered normalizer |
| constant h=f、off Z=N | 指定既約D源は空 | 全K1 bad u-row零、actual u因子 |
| constant h=L^3 | 既essential-threeの八23成分へ全被覆 | Ore全18必要源、一つの元Jacobi式 |
| moving degree一、全t nonzero | 像23以下 | QQ(t)及び全七exceptional walls |

全leading還元は `PLANE_M2_DOUBLE_ROOT_CONNECTION_GEOMETRY.md`、
`PLANE_M2_MOVING_DOUBLE_ROOT_NORMAL_FORM.md`、
`PLANE_M2_DOUBLE_ROOT_CONIC_REDUCTION.md`、
`PLANE_M2_DOUBLE_ROOT_CONSTANT_NORMAL_FORMS.md` 第2、3節。
primitive rootの非primitive境界はconstant rootへ戻し、t=0はnormal位数三以上へ明示する。

有限source証明は `PLANE_M2_DOUBLE_TANGENT_FINITE_BOUND.md`、
独立別算術/normalizer査読は `PLANE_M2_DOUBLE_TANGENT_SEPARATE_REVIEW.md`。
producer reportのwhole-Gフラグは粗上界27が未閉だった当時のまま保存する。
このフラグを修正して証明を代用せず、次の構造補題で残源を別に閉じる。

有効なOre証明は **`PLANE_M2_DOUBLE_TANGENT_CUBE_ORE_RECOGNITION_V2.md`**。
全K1九及びP2全九でB=C[z,u,v]が可換、{y,B}=delta(B)、{x,B}=yA(B)+V(B)。
元quad(x,z,u,v)はactualDのy微分の符号反対なので、全solutionでactualDはBだけに依存する。
partial_x actualD零も全sourceで恒等的。
既約rank四なら既essential-three八23成分定理が全pointへ適用出来る。
独立別実装の全18source検証は `PLANE_M2_DOUBLE_TANGENT_CUBE_SEPARATE_REVIEW.md`。
未受理V1のstrict-Jacobi誤段は不採用で、ここへ依存しない。

この定理はGのwhole非優越であり、Z0及びmoving全sourceの個別E帰属を追加しない。
normal-doubleのtransverse conic+lineは別独立証明で全空。
triangleの四leading源、concurrentの三leading源のfullhigher、nonreduced R及びm3は別残余。
全m2 squarefree reducedRの既全MPS被覆も本定理とは別に保持する。
