# D5+tangent 全source：唯一の正当なanchorでの全周囲結果

2026-11-07。PRIMARY E。全marked family・全六磁気源の還元は
`PLANE_HYPERPLANE_QUARTIC_D5_TANGENT_SOURCE.md` による。
そのcanonical R=1点一つについて、許可された最弱の全周囲計算を行った。

`classification_ga9_cuspidal_quartic.py --d5-source` は既存runnerのisolated modeであり、
元の全FunctionDef ASTを旧公開producer bytesと照合して全て不変であることを確認する。
旧main/defaultを再実行・変更していない。

全39/90の部分sliceではなく、全150 coefficient・全350 projective Jacobiから

```
rank D=101, Euler vertical=24,
rank full51 Euler/GL/scale=48,
rank full51+six source derivatives=48,
D*all_actual_directions=0
```

をQQで直接確認し、各非零小行列を保存した。
全六磁気源のstrict10 Jacobiも恒等的に零で、実五次は16hF。
よってsourceの実23軌道下界を再確認し、全周囲の接空間上界は射影24である。

独立な48 symmetry directionsを消す真の102次元linear cutを構成した。
そのkernelは正確にnormal一本であり、48+1+101が全150空間の直和となる非零小行列を保存した。
唯一normalの全350 quadratic Jacobiは非零だが、

```
rank[D | quadratic_Jacobi(normal)]=101=rank D
```

であり、**二次余核は零**である。
これは原点根基・成分性・実曲線の存在を証明しない。
契約に従いこの段階で止め、cubic・higher IFT・別anchorを計算していない。

結果は `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_D5_AMBIENT_20261107.json`。
finite job `jobs/plane_hyperplane_quartic_d5_ambient.json` はgate PASS、
120秒/1536MiB契約、実測5.139秒/65,216KiB。
全行列とnormal/quadratic vectorを既存basis順で保存している。
独立arithmetics review及び全既知成分との相違監査はまだ未実施。

現時点で証明できるのは全D5+tangent源が一つの23軌道閉包であることと、
この点の接空間dimension24までである。全Poisson成分か既知の大きい成分の境界かは不明。
同じnormalの次数を上げる次手は採らず、実幾何族又は別の全marked branchを次の対象とする。
