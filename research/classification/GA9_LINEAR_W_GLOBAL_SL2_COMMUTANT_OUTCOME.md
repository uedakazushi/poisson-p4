# Ga9 global SL2＋全九不変補正族：一回の全係数unit障害

2026-11-08、PRIMARY E、**ONE ansatzのexact陰性**。T37全非所属、Ga9新成分性又は全E帰属は主張しない。

## 1. 全源と一回だけの実行

`GA9_LINEAR_W_GLOBAL_SL2_COMMUTANT_FAMILY.md` の全projective invariant basis九本と、真のcubic product係数r0,…,r3を保持した。

```
qε=(s+εt)(s−εt),
cε=Σ r_j(s+εt)^(3−j)(s−εt)^j, Σr_j=1,
Πε=ρ2(qε)+ρ3(cε)+εΣ_(i=0)^8 k_i K_i.
```

rootのpure basis／dim13−4／module normalization／全period査読と明示GOの後、別の実行caller・jobを新identityで凍結した。準備契約db99acea…のbytesは変更しない。

14入力hashと180秒／3GiB gateを確認後、一回だけ実行した。全150係数のpolynomial族、全350 Jacobiのε次数0–6、全五Pfaffian及びsigned quintic全126係数を構成・保存した。全SL2 invarianceと有限module恒等式、endpointも元係数で検証した。

実測1.161948秒、maximum RSS54332KiB。新point、ambient operator、rank／minor探索、normal profile又は高次jetは無し。Groebner basisを計算する必要もなかった。

## 2. 全13未知定数からの可逆線形消去

全元Jacobi係数とΣr−1の線形部分のrankは九。元行への完全なrow-operation liftを保存し、次の全解を得た。

```
k0=−9k8/4, k5=−9k8/2,
k1=k2=k3=k4=k6=k7=0,
r0=1−r1−r2−r3.
```

従って残未知数はk=k8,r1,r2,r3である。最初からこの一dimensional Kへsourceを制限したのではなく、全九kと全四rから消去した。

Kは定数倍を除いて

```
D1∧D2+2R+∧R−−(4/9)w²J3
```

へ制限される。これは前の純手centralizer計算とも一致するが、算術はその手縮約を仮定していない。

## 3. 四つの元係数だけで生じるunit

保存350-row順序でE_(i,d)をJacobi row iのε^d係数とする。全線形解を代入した元式のうち、以下はそのままの係数である。

```
E_(100,2)=32kr1+64kr2+96kr3−48k+(256/3)(r1+r2)+192,
E_(106,4)=−36kr1−36kr3+18k−(32/3)(r1+r2)+8,
E_(115,2)=−36kr1−72kr2−108kr3+54k−(64/3)(r1+r2)+24,
E_(118,4)=288kr1+288kr3−144k−(256/3)(r1+r2)+64.
```

いずれもquadruple(h,x,y,w)の式であり、coordinate monomialは順にhyw²、x³y、xy³、xyz²。
定数乗数だけの恒等式

```
(9/2368)E_(100,2)+(7/592)E_(106,4)
 +(1/296)E_(115,2)+(7/4736)E_(118,4)=1
```

を全四残未知数のpolynomialとして厳密に検証した。
例えば9E_(100,2)+8E_(115,2)は `1792(r1+r2)/3+1920`、
8E_(106,4)+E_(118,4)は `−512(r1+r2)/3+128` となる。
第一式と第二式の7/2倍を加えると定数2368である。

したがって全normalization／全350係数の解は存在しない。有限体、標本の残未知値、selected jet又はdegree増加で得た陰性ではない。
線形消去の完全な正逆と元行liftがあるので、quotient内のunitは元の全13定数sourceの不整合を証明する。

## 4. exact evidenceと射程

report `CLASSIFICATION_GA9_LINEAR_W_GLOBAL_SL2_COMMUTANT_20261108.json` のSHA256は
`0007109f8a0277d50321523160b82d225cd0108e5191c00521d9430ed31dc2c1`。
full source／全2450 Jacobi coefficient slots／全Pfaffian／全linear quotient式／unit lift／inline sourceを保持する。
実行source SHAは `852a6a156e6c6ed6c77d2367b8debdf978994dd0afd64dfef30a1d426ee790d7`。

今回閉じたのは**定義したεKとmoving-eigenframe cubicからなる、この全九補正ansatz**だけ。
全T37 arcへの必要次数定理を仮定していないため、他のglobal invariant構造／normalization residue族まで排除しない。
ε²K、追加module、別normal、別point又は高次jetでこのansatzを繰り返し延長しない。

fixed-D全源とunique singular plane rank二によるGa9認識定理は有効のままである。
次のprimary義務は別のglobal Poisson-period extension又はwhole parentへ効く幾何であり、Ga9 E帰属は未解決。
catalogue、主稿、global台帳とGitはこの作業で変更していない。
