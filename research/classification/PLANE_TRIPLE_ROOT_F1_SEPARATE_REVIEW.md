# m4 generic31：全 F1 必要源は共通 vertex 外積、rank四空

2026-11-07、別作者 `/root/quartic_review`。**PASS**。
root の新全 L4 polarizations を、別の Fraction・疎 Laurent 多項式実装で再構成した。
producer/CAS は import しない。新 anchor、jet、GB は無い。
有限契約 `jobs/plane_triple_root_F1_separate_review.json` は実行前 gate PASS、
全1075×150 再構成を含め 0.982秒、39,460 KiB。

## 全元の方程式から完全 L4 へ

座標は (A,B,u,v,w)、実 quintic は

```
F1=u4 A+u³v B+u²v³+v5.
```

全150 quadratic skew係数に対し、`Π#dF1=0` の5×210係数と divΠ=0 の25係数を
ゼロから構成し、旧 report の全1075×150行列と一致させた。
保存 rank146 minorを Fraction で独立再計算し、全四 kernel の独立性及び全1075行との積零を検証した。
従って必要源 L4 は完全であり、chosen slice や少数方程式の kernel ではない。

この全 kernel は実は次の四本 `Zi∧∂w` そのものである。
Zi の成分順は (A,B,u,v)。

```
Z0=(-uv, u², 0, 0)
Z1=(-v², uv, 0, 0)
Z2=(-4Au-5Bv, 2Bu-2v², u², 0)
Z3=(-5Av-Bu-3v², 15Au+(15/2)Bv, -(5/2)uv, u²).
```

全係数 c0,…,c3 を保持した任意の源は

```
Π=(Σ ci Zi)∧∂w,   ∂w Zi=0.
```

したがって全五 raw Pfaffian は**恒等零**、全 strict Jacobi も構造的に恒等零である。
独立実装では更に全十 quadratic parameter monomial に対し、全350 projective Jacobi行と
五 Pfaffian の全70四次係数を直接展開し、root の保存全矩形と一致させた。
全 actualD 五次係数も零となる。全四係数が零となる境界を含め、局所化は一つもしていない。
従って **F1を actual rank4退化因子に持つ Poisson源は空**。
全Jacobiが零であることを rank4非空性と誤読しない。全 source は rank≤2である。

## 全 generic31 actual D への正当な橋

`PLANE_M4_REPEATED_SEPARATE_REVIEW.md` の全幾何を再読した。
image-rank二の generic31 leading は fixed triple root の Bertini 補題により

```
D=u³(Au+Bv)+g5(u,v),  A,B independent.
```

既約性により κ=g5(0,1)≠0。κ0は D がuで割れる別枝で、既約Dの対象外である。
A,Bへ normal一次式を加える可逆 linear shear により u5,u4v,u3v²を吸収する。
残る `αu²v³+βuv4+κv5` の βは `v→v−βu/(5κ)` で消し、
再度 A,B shearで生じた高u項を吸収する。
新係数は `α−2β²/(5κ)`。κの非零尺度、uの有限平方根尺度により、
残係数零/非零の**両補集合**を F0/F1 へ戻せる。
これは true GL5 変換であり、leading connection、Rの被約性、logpoleを仮定しない。
両 F は Aに関する primitive linear polynomial なので絶対既約で、essential-four。

canonical divergence-zero strict quadratic lift は、全 projective tensorから一意に取れ、GL5-equivariant。
BV恒等式及び Cartan による Pfaffian–Casimir `Π#dD=0` は essential-three に限定されない。
実際 divΠ=0 の Poisson liftでは `ι_(Π²)Ω` が閉一次形式となり、homogeneous Eulerから定数倍dD。
rank4域で Π がこれを殺し、恒等式として全域へ延長する。
従って actualD class[F1]を持つ全 projective rank4源は、上の完全 L4 に入る必要がある。
Euler gauge をもう一度 dimensionとして引いたり、任意代表を div0と仮定していない。

F0の全源は既別作者証明で rank≤2。
その全 quadratic gradient syzygy は S1Z+S1Y、AR1=span(Z,Y)、AR0=0で、
vertex方向を戻しても skew imageは最大三次元、従って偶数rank≤2となる。
F1の今回の空性と合わせ、**R非被約も含む全 generic31 image-rank二 rank4源は空**である。
旧F1 G像≤21上界は偽にはならないが、必要支持としては今回の空性が強い。

image-rank一/generic22/4は別の essential-three八源定理へ送る。
211、squarefree、m2/m3及び全Eの結論をF1のこの一点正常形計算だけから推論しない。
他の independently accepted m4全枝との合成は integrator の別 aggregate命題で行う。

全入力と実行 reviewer bytes は
`CLASSIFICATION_PLANE_TRIPLE_ROOT_F1_SEPARATE_REVIEW_20261107.json` と対応manifestに固定した。
旧 producer report、旧 G≤21証明書、root report は変更していない。
