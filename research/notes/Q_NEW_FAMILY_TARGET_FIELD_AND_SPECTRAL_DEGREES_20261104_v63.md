# 第一jet像 B_NEW の不変関数体の指数1と spectral degree 10・20・40

2026-11-04。これは既知の有限 source から **実際の第一jet target** への関数体包含を
尽くすための、独立に監査できる帰結である。ZN 原始元の独立 gate が未受理なら
本定理は条件付きであり、専用報告の `status` と `pending` を参照する。
受理後に限り指数1と下記の次数を無条件の帰結として採用する。
巨大な norm、Inv14 の新展開、source26、新しい点・prime・CAS は使わない。

## 1. 実際の関数体と入力 gate

まず係数体を `C` とする。既受理の既約15次元 primitive first-jet hypersurface
`B_NEW subset W16` の Sp4 不変有理関数体を `E=C(B_NEW)^Sp4` と書く。
**B_NEW は25次元 quadratic Poisson componentそのものではない。**
既存の実 source から B_NEW への支配写像は、
全 Sp4 invariant の pullback を source derivation の定数体、さらに既知の有限
involution の固定体へ送る。その**実際の写像に束縛した包含**を

```
E subset L = C(xi,r,s,q,B)
```

と使う。任意の hypersurface の invariant fractions に関する一般論から
この包含を推論するのではない。

既受理の両向き Möbius 回収により

```
L=C(xi,R,Z0,r,s),       Kshape=C(R,Z0,r,s),
K0=C(R,Z0,Z1,Z2),       L=Kshape(xi),
```

であり、`xi` は shape 四変数から独立である。全て表示された非零分母の一般部に
おける関数体の主張であり、除外 boundary の延長・properness は不要である。

実際の polynomial invariant を `F4,p0,...,p4,N12` と記す。degree はそれぞれ
4,8,...,8,12 で、元の規格化のまま

```
F4²=400p0+1600p4,
R=(p3+4p4)/p0,
Zi=pi/F4²  (i=0,1,2),       ZN=N12/F4³.
```

従って `K0` と `ZN` は `E` に入る。octic base への実際の優越性から
`p0,...,p4` は代数的独立であり、`R,Z0,Z1,Z2` も独立である。

必要な独立 gate は次の通りで、専用 receipt が全 report hash を固定する。

- SOURCE_CONSTANT_FIELD と LOW_INVARIANT_BIRATIONAL_COMPRESSION:
  上の実 source 包含と既知固定体の五つの有理座標。
- TWO_SHAPE_RECOVERY: q,B の両向き有理逆と正確な normalized invariant maps。
- OCTIC_BASE_DOMINANCE: 五つの p の独立性と規格化。
- TWO_SHAPE_CONTENT_PRIMITIVITY と COPRIME_DEGREES:
  `Kshape=K0(r)`、`[Kshape:K0]=10`、s の有理回収。
- GENERAL_SOURCE_SCALE: 特定の低不変式からの外挿でない、全16座標による
  任意次数の cocharacter covariance。
- FULL_DEGREE24_SOURCE_RELATION: B_NEW の斉次 prime equation の degree は少なくとも26。
- ZN_PRIMITIVE_KRYLOV の **独立レビュー**:
  `K0(ZN)=Kshape`。有限自由 rank10 係数行列の非零 determinant 21349
  modulo32003 が証明する範囲であり、有限体の値だけを QQ/C の消去結果に昇格しない。

## 2. 次数10は C 上でも正確である

ここで QQ 上の既約性だけを C 上へ無条件に移してはならない。
既存の証明はその不足を持たない。`C(r,R,Z0)` 上でも、s の二つの有理写像の
次数は3と5で、保存した numerator/denominator の coprimeness 恒等式は保たれる。
よって joint map は birational、primitive implicit equation はその既約曲線の一乗である。
content の最大性に使う rational-coefficient Bézout 恒等式も C(R) 上でそのまま有効。
従って全 content 除去後の G は C 上でも幾何的に同じ incidence の prime equation である。
Gauss lemma と r-degree10から

```
[Kshape:K0]=10.
```

ZN の Krylov determinant は、monic G による **有限自由 rank10** algebra の
係数 determinant の非零性を示す。この determinant が QQ で非零なら C でも非零。
generic denominator が非零という受理済み条件の下で

```
Kshape=K0(ZN) subset E.                           (1)
```

特殊 algebra で分母の逆元を仮定する必要はない。これが target への包含に必要な
shape 全体の回収である。

## 3. 奇数 xi-weight の不変式による xi の回収

受理済み全16座標 cocharacter は、source の c を lambda² 倍する操作が
beta の lambda 倍と Sp4 作用との合成であることを示す。`xi=40c/(Bu)` であり、
shape 座標変更は c に依存しない。従って任意の斉次 degree2n invariant は

```
I_(2n)=xi^n h_(2n)(shape).
```

特に `F4=xi² f4`、`f4` は Kshape の非零元である。
非零の普遍 degree14 invariant `J` を一つ選ぶことができる
（既受理の Inv14 次元は3）。B_NEW 上で J が零なら、その高さ1の斉次 prime
equation が J を割る。degree 下界26と `degree J=14` に反するので、J の制限も非零。
従って

```
J=xi^7 h14,       0!=h14 in Kshape.
```

(1) と F4,J が実際の invariant であることから

```
xi = (J/F4³) * f4³/h14 in E.                     (2)
```

既に E subset L なので、(1)–(2) は

```
E=L=C(xi,R,Z0,r,s),        [L:E]=1.              (3)
```

を証明する。特に **B_NEW の実際の orbit-invariant field は純 C 有理・超越次数5**。
J,h14 の巨大な式を展開しなくても、存在・非零性と正確な field inclusion で十分である。
これは明示的な expanded inverse formula を保存したという主張ではない。

## 4. F4 を残した affine spectral base への正確な次数20

`H=C(F4,R,Z0,Z1,Z2)=K0(F4)` と置く。`F4=xi² f4` で xi は shape から独立なので、
F4 も Kshape 上 transcendental。従って

```
[Kshape(F4):K0(F4)]=10.
```

`xi²=F4/f4` の右辺は Kshape(F4) の F4-adic valuation が1であり、平方ではない。
従って `T²−F4/f4` は既約で、標数0なので分離的な二次拡大となる。ゆえに

```
[E:H]=[L:Kshape(F4)]*[Kshape(F4):H]=2*10=20.      (4)
```

ここで任意の二次 base change との linear disjointness を仮定していない。
transcendental F4 を先に加えることと、valuation による二次既約性を別々に使う。

## 5. 五つの octic coefficient base への正確な次数40

`P=C(p0,p1,p2,p3,p4)` と置く。次の有理逆式により `H=P(F4)` である。

```
p0=F4² Z0, p1=F4² Z1, p2=F4² Z2,
p4=F4²(1/1600−Z0/4),
p3=F4²(Z0(R+1)−1/400).
```

P は独立な五変数の有理関数体。既約な linear polynomial `400p0+1600p4` の
valuation は1なので平方ではなく、`F4²=400p0+1600p4` は正確な二次拡大を与える。
(4) と塔の公式から

```
[E:P]=20*2=40.                                  (5)
```

これは generic affine invariant quotient の field degree。全 base 上の有限・flat
morphism、未調査 boundary の fiber 数、ramification の全分類は主張しない。

## 6. Projective normalized spectral degree は10

beta の scalar lambda 倍は、Sp4 を除いた既知 source 座標で
`xi -> lambda² xi`、shape 固定として作用する。C 上で lambda² は C* 全体を走るので、
(3) の scalar-fixed field は正確に Kshape。従って projective B_NEW orbit quotient の体は

```
E_projective=Kshape,
[E_projective:K0]=10.                            (6)
```

ここで「projective」は beta の scalar 同値をさらに取る意味である。
既知 source の有限 involution をもう一度二重に数えたり、affine xi の二次拡大を
projective field に残したりしない。

## 7. QQ 上での定義体監査：同じ指数1と有理性

上の結論は、この B_NEW に関しては QQ まで降りる。これは ambient W16 商の
Shepherd-Barron の C-rationality を QQ-rationality と読み替える推論ではない。

元の六変数 source、全 primitive beta 有理式、Sp4 表現、source kernel、section、
既知 tau の実際の有限同値は全て QQ 上に定義される。source の image closure
`B_NEW,Q` は QQ 上の有理 source の像である。kernel は係数体の flat extension と
両立し、C 上の同じ像は既約であるから、この Q-image は幾何的に既約である。
従って

```
E_Q=QQ(B_NEW,Q)^Sp4 subset L_Q=QQ(xi,R,Z0,r,s)
```

という実際の pullback inclusion がある。derivation の定数体計算も既知有限
involution の固定体計算も、QQbar の点や非有理な group element を必要としない。
後者の実 group matrix は QQ 上の diagonal sign matrix である。

G、content の全 Bézout、linear-C1 回収、ZN の有限自由 determinant は全て
QQ 係数のため、`Kshape_Q=K0_Q(ZN) subset E_Q` がまず QQ 上で成立する。
Inv14 の定義条件は rational Lie matrices による有限な QQ 線形方程式である。
線形 kernel は C への拡大と両立するので、既受理 `dim_C Inv14=3` から非零の
QQ 係数 invariant J14 を選べる。C 上の degree 下界26により、その QQ-source
制限も非零である。全16座標 cocharacter の有理恒等式は QQ 上でも成立し、

```
J14=xi^7 h14,     f4,h14 in Kshape_Q,     f4*h14!=0.
```

従って (2) の xi 回収は QQ の有理演算だけで成り立ち、

```
E_Q=L_Q=QQ(xi,R,Z0,r,s)
```

が得られる。**B_NEW の orbit-invariant field は純 QQ 有理・超越次数5**である。
QQ 上の spectral field degree も同じ10・20・40。二次既約性は同じ奇数 valuation
で示され、projective scalar-fixed field は weight2 の代数的 Gm 作用の固定体である。
QQ の平方が QQ* を全射するという誤った主張を使わない。

この QQ 結論は、ambient `QQ(W16)^Sp4`、25次元 Poisson component 全体、
その source-to-component field index、あるいは全射的な ambient invariant fractions
に関する結論ではない。

## 8. 範囲と保存

本ノートが閉じるのは **第一jet像 B_NEW の actual invariant target と既知 fixed source の指数**および
normalized projective/二種類 affine spectral base への正確な generic degree である。
ambient W16 の全 invariant field に関する既存の rationality theoremとは別の主張である。
curve genus、actual global Poisson component の rationality、NEW hypersurface の
具体的な定義 polynomial・その正確な degree、全 Poisson component 数・exhaustion は
ここでは決定しない。いずれもこの field-degree 計算と混同しない。

専用 checker は入力レビューと hash を結ぶだけで、新しい数式計算は行わない。
独立 ZN gate が無い場合は通常実行を停止し、`--allow-pending` の明示指定時のみ
結論を `null` にした conditional report を保存する。最終 receipt は全 gate 通過後に
新しい output へ保存する。

```
.venv/bin/python code/q_new_family_target_field_and_spectral_degrees_20261104_v63.py \
  --output /tmp/new25-target-field-degree-replay
```

既存の凍結済み「指数未決」ノートは、その時点の正確な記録として変更しない。
