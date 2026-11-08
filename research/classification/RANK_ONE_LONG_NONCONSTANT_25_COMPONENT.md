# rank一long・非定数gの全源閉包：新しい25次元成分

2026-11-08、PRIMARY E_multiple_hyperplane。**producerの全源成分証明。独立fullambient reviewは別receipt待ち。** 「新しい」は固定prior129に対する相違を意味し、外部文献全体への新規性確認と同一視しない。

## 1. 全源の下界はpointに依存しない

`RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md` の全源Sは非定数affine g、全affine U with U(g)=(div U)g、全二次Kを保持する。係数から全parameterを一意に復元でき、Sは既約12次元。正確なAff4 modular stabilizer十次元、scalar stability及び唯一double-H incidenceにより、そのPGL飽和射影像閉包Zは**正確に25次元**である。

同noteのfullprior129比較は、既知成分の全rank二閉条件、dimension、proper factor incidence及び三つのMPS25の被約開によりZがそのいずれにも含まれないことを証明する。入口の独立純review `RANK_ONE_LONG_NONCONSTANT_ENTRANCE_SEPARATE_REVIEW.md` はこの全source／全fiber／全priorの橋を受理した。以下の一点からsource次元又は全帰属を推測するのではない。

## 2. 唯一実行した固定点と全空間

固定X=∂x+x∂yの元sourceから、次の一つの点だけを使った。

```
g=z,
U=2z∂z+(3+5z)∂w,
K=1+2z+3w+5z²+7zw+11w²,
F=(U+3)K=12+46z+75w+70z²+145zw+33w².
```

原点の評価をsymplecticに固定する操作や追加sliceは課していない。h-Casimir liftを二倍して整数化した。全十strict Jacobi／全350 projective Jacobiは零、modular場は2h∂x+2x∂y、実退化因子は

```
D=−4h²*z*Q,
Q=h*y−x²/2+h²+2hz+3hw+5z²+7zw+11w².
```

Qのsymmetric matrix determinantは非零で、rank五である。全150二次lift係数における全350式の微分D、全Euler25列V、全Euler／GL／scalar51列Mを構成し、QQ上で

```
rank D=100, rank V=24, rank M=47, D*M=0
```

を得た。全matrix、選択nonzero minorとその値、全kernelとfree-coordinate identity、全五つの4×4 Pfaffian、全signed実Dを保存している。rank Mから得る**この点**のPGL軌道次元は22であり、それを全source又は全成分の次元とはしない。

有限契約は `jobs/rank_one_long_nonconstant_one_ambient.json`。gate PASS後のONE runは2.098528秒、最大RSS55,744KiB。結果は `CLASSIFICATION_RANK_ONE_LONG_NONCONSTANT_ONE_AMBIENT_20261108.json`。同じ点のselected sourceだけでなく、全350×150 ambientを使用した。別点、追加profile、jetは実行していない。

## 3. 下界25と全接上界25の一致

Euler商前の150係数空間では、Zのaffine coneと全Euler24の逆像は既約50次元であり、この実点を含む。従って全Jacobiスキームの局所次元は少なくとも50。一方rank D=100なので全Zariski接次元は150−100=50以下。局所次元と接次元が一致するから局所環はregular、唯一の全空間成分を通る。

Euler24とscalar一を除くと射影局所次元は25。閉既約Zも25次元なので、この全成分に等しい。これで**全long非定数g sourceの閉包Zは新しい25次元全Poisson成分**となる。sourceの各pointのsmooth性、局所唯一性、generic PGL orbit又は全nonreduced構造を主張していない。

## 4. 全境界とshort源への橋

非定数gを保った全Sは同一の閉成分Zに入り、UやKの全零係数／全resonance／quadric rank低下も含まれる。選んだ点の非零係数を全Sへ条件として戻さない。

`RANK_ONE_SHORT_TO_LONG_WHOLE_POLYNOMIAL_CURVE.md` の分母のないactual二次Poisson曲線により、**全short非定数g源もZに入る**。λ非零では真のaffine変更で上の全long源へ移り、λ零で元shortの全六係数に戻るためである。非線形rectificationをPGLと取り違えた推論ではない。

同curveは全short定数g／H3零を既存NEW25へ送る。一方、long g零全Fとlong定数gの三次共鳴は未解決であり、それぞれのshort補集合も対応するlong残余へ戻した段階である。これらをZへ含めた、又はrank一c非imageの全支持が完成したとは主張しない。
