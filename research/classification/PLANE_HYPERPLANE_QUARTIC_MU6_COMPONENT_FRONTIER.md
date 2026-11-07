# μ6/r2全十一源：九component候補、旧DPMへの全源認識、一残余

2026-11-07、PRIMARY E。producer fullambient/cut証明、独立確認前。
元の全十一marked source・全六磁気の一23軌道閉包還元は別二noteによる。
本結果は新しい九成分を即時計数するものではない。

## 全十一の一回のambient batch

全150係数・全350Jacobi、Euler/GL/scale全51方向と全六source微分をQQで構成した。
各canonical Rはcokernelの最初の非零係数で選んだ一代表のみ。
全十一でEuler rank24、作用とsource rank48、D*全実方向=0。
overall primitive integer scaleはprojective同じ点で、別anchorを選んでいない。

| source | rank D | 最弱の全周囲結論 |
|---|---:|---|
| A6 cusp-secant/bitangent |102| smooth unique23 component候補 |
| A6 smooth flex |101| 真一本normalの非零二次余核、23component候補 |
| A4+A2 cusps join |102| smooth unique23 component候補 |
| A4+A2 A2 tangent |101| 非零二次余核、23component候補 |
| A4+A2 A2-secant/smooth-tangent |102| smooth unique23 component候補 |
| A4+A2 smooth flex |101| 非零二次余核、23component候補 |
| 3A2 cusps join |102| smooth unique23 component候補 |
| 3A2 cusp tangent |101| 二次余核零、局所反復停止 |
| 3A2 smooth bitangent |102| smooth unique23 component候補 |
| E6 smooth flex |100| 二normal、局所反復停止 |
| E6 smooth bitangent |101| 非零二次余核、23component候補 |

D102ではraw tangent48が実orbit48と一致し、smooth projective23/unique componentとなる。
D101の四positive caseでは、実symmetry48に対する真の102次元cutを取り、
48+normal1+implicit101が全150の直和である証明書を保存した。
全350二次normalベクトルに対し、全150 D列を消すfunctionalの値が非零。
よってIFT後の唯一normalの残式は非零二次leadingを持ち、被約cutは原点。
Krull上界raw48が実orbit48と一致し、全23軌道閉包が全空間componentとなる。
非被約厚さや高次jetは不要で、smoothnessはD102だけに主張する。

report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_AMBIENT_20261107.json` はCOMPLETE。
59.501秒、最大RSS79,896KiB。D5+tangentは再計算していない。
二nonclosing caseはcubic/higher IFT・別anchorへ進めず、そのまま保存した。

## 旧actual witnessとの正確な比較

旧五つの23/H+quartic成分をmetadata名だけで処理しない。

* ESSENTIAL5：`legacy/v06/nilpotent22.tex` の実五次はhF、
  `F=h4+h²xz+x³y−x³z+x²z²/2−x²zw`。
  Dの方向微分のrows `(h4,h³z,hx³,h³x,hx²z)` は、columns `(h,x,y,z,w)` でdet−5のminorを持つ。
  従って実essential rank五である。新source全体はclosed essential≤3条件に入るため等しくない。
* MINIMAL_M：印刷された六括弧から実Pfは`h4+3h²x²+12x³z`、発散は`6x partial_y`。
  旧23orbitはrank一canonical cofactor、新十一は全rank二なので実orbitは異なる。
  旧componentはその稠密23orbit閉包と受理済みである。
* A0：印刷六括弧の実Pfは`h4+h²xz−z4/4`。
  Fは唯一D5、marked H=hはそのcusp tangentである。
  これは新mu六候補の全cuspidal profileと異なる。
* DPM：`legacy/v07/dplus_fullp4.tex` の実Fは`h4+h²xz+p_D(x,z)`。
  `p_D` のroot multiplicityは(3,1)、Fには三A2、H=hはその一つのcusp tangent。
  新九候補の3A2はHが二cuspを結ぶ又はsmooth bitangentであり、marked geometryが異なる。
* Ga9constant：受理済み全十sourceのactual canonical cofactorは`y partial_x+h partial_z`、
  marked H=hはdual cofactor image内。
  新全sourceではhはimage外である。
  同じ23componentなら、旧dense constructible actual source imageと新dense23orbitが交わるので矛盾する。
  image-membershipのclosed性を無証明で仮定する議論ではない。

旧actual六括弧と上の発散は `GA9_CONSTANT_ANCHOR_REVIEW.md` と対応独立checkerが別入力で検証済み。
essential5/A0のPf式も印刷sourceから直接得られ、cofactorの名前から型を仮定していない。

## nonclosing sourceの正の旧成分認識

DPMの旧点は全51作用rank48、実PGL orbit23であり、旧componentはそのorbit閉包。
quarticの三cuspとmarkedH=cusp tangentというpairは、今回の全tricuspidal normalizationにより一orbit。
旧六括弧はbase-bracket零のcotangent sourceで、anchor determinantもそのF。
全quadratic Saito spaceは二次元なのでbase変換後の旧anchorは新basisの定数GL2結合。
全六磁気も可逆pullbackで移る。旧orbit23ならmagnetic invariantは非零でなければならない。
零ならfiber homothety stabilizerが一つ増え、全orbit23に反するためである。
全mag normalizationにより旧点は新canonical actual orbitそのものに入り、
**3A2 cusp-tangent 全六source閉包は旧DPM componentに等しい。**
二次余核零を理由に旧cubicを反復する必要はない。

同じ全pair/Saito/磁気議論で、既済D5+tangent canonical sourceは旧A0の実orbitと等しい。
旧A0のPfは上記D5 pair、旧orbit23によりmagnetic invariant非零。
従ってその全源も**旧A0**へ被覆できる。
D5のD101/Qzeroを新componentや未帰属とする必要はない。
これらはcofactor型だけからmembershipを推論したものではなく、全marked pairの可逆正規化と
完全なSaito anchor/magnetic sourceを経由したactual orbit一致である。

## 全prior111及び相互の相違

priorは公開commit37d56e2a、catalogue blobdfc7c037b2050012fbe65c6c34b0e1ad53d3ff59。
完全111行を `inputs/mu6_prior111_37d56e2a_20261107.json` に固定した。
rank二又はdimension23以外とはgeneric rank/全component dimensionで相違する。
同23のMPS29は受理済みactual toric-SNC点を持ち、実Dは非零κ*hxyzwなのでessential rank五。
新閉包全体のessential≤3というclosed条件から相違する。
残る旧非MPSの異factor degree・nonreduced・irreduciblequinticはproper multiplication locusで処理する。
同じreduced1+4の旧五つだけは、上の実witness比較を使う。

新九候補相互は全sourceがdense23orbitなので、同componentなら実orbitが交わる。
quartic singularity profileが違う場合は不可能。
同profileではA6のHがcuspを通る/避ける、A4+A2のH上のcusp数・contact、
3A2のH上のcusp数二/零で区別できる。E6positive候補は一つだけ。
従って独立にcomponenthoodが受理されれば九つは相互及び全prior111と異なる。
まだrootのglobal計数は変更せず、別作者の全geometry/arithmetic/audit reviewへ渡す。

E6 smooth-flex sourceは接上界25、実source23まで。旧25/24親又は新componentかは不明で、
actual known-family認識など別の幾何経路を要する。二normalの次数を上げない。
