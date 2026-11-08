# unimodular対数全境界：零一次被覆と一つの残実軌道

2026-11-08、PRIMARY `E_multiple_hyperplane`。producer純証明、独立/root査読前。
旧Gの全kernel/minor/全residue証明書を再計算せず、その**全残Pf式**からE membershipを導く。
新しいoperator・点・minor・jet・GB・CASのrunは無い。

## 1. 入力と零一次既知閉和

旧 `UNIMODULAR_LOGARITHMIC_SUPPORT.md` と
`UNIMODULAR_LOGARITHMIC_POLE_BOUNDARY_REVIEW.md` は、全最高部と全C6/L24について
translation image、真のnormal slice、六定数の可逆消去、全残Pf式と全minor分母の根を保存する。
全translationの逆を最後まで保持する。Pf nonzero constantなら任意のtranslated基点でCも非退化。

L=0を得たP=C+Qはactual linear Darboux frameでJ+Qへ移る。
元の全Poisson条件は[J,Q]=[Q,Q]=0なので、
`BETA_ZERO_COVERAGE.md` のQ45全被約22源から既知閉和への既証明を適用出来る。
その有限閉PGL和を以下\(\mathcal B\)と書く。
Q22自体を本稿の数式整理で再証明したとはしない。

## 2. 実L112・full mixed Mの全residue

最高部は順序(xy,xz,xw,yz,yw,zw)で

```
Q_u=(0,−xy,x²,u y²,−u xy,(u−2)xz+(2u−1)yw).
```

旧全source証明により、全u（generic u=2/3のforcing補集合を含む）でn0=n1=0。
u=0、1、infinityは後述の既別型へ戻す。
既saved full残式にn0=n1=0を代入すると、generic uの全normal sourceは

```
L=(0,0,n2*x/(u−1),n3*y,−u*n2*y/(u−1)−n3*x,n2*z+n3*w),
C13=−R/[u(u−1)], C24=a0,C34=a1, C12=C14=C23=0,
R=a0(u−1)+n2*n3,
Pf(P)=a0*R/[u(u−1)].                              (1)
```

全残式のうち二つは−n2*R/(u−1)²、−n3*R/[u(u−1)]。
Pf非零からa0R非零なのでn2=n3=0、従って全L零。
この結論はu∉{0,1,1/2,3/2,2}上の全QQ(u) sourceの恒等式であり、数値period標本ではない。

例外residueの全saved kernelでも同じunit結論を直接使える：

| u | Pf nonzero constant | 全残式に含まれるunit-forcing | 結論 |
|---|---|---|---|
| 3/2 | 2a0(a0+2n2n3)/3 | n2(a0+2n2n3), n3(a0+2n2n3) | 全normal L零 |
| 1/2 | 2a0(a0−2n3n4) | n3(a0−2n3n4), n4(a0−2n3n4) | L=n2*y∂z∧∂w |
| 2 | a0(a0+n3n4)/2 | n3(a0+n3n4), n4(a0+n3n4) | L=n2*x∂z∧∂w |

すべてn0=n1零後の全sourceであり、追加normal変数を捨てていない。
最後の二行も次のliteral full familyで\(\mathcal B\)へ入る。

```
C_u13=−a0/u, C_u24=a0,C_u34=a1, 他0,
P_u^y=Q_u+C_u+n2*y∂z∧∂w,
P_u^x=Q_u+C_u+n2*x∂z∧∂w.                         (2)
```

a0非零、u非零上のactual quadratic familyである。
generic uでP_u^yはQ_u+C_uのargument wをw+n2/(2u−1)へ、
P_u^xはargument zをz+n2/(u−2)へtranslationしたもの。
Q_uの他五括弧はz/wに依存しないので、これは全六括弧のliteral等式である。
従ってgeneric uで\(\mathcal B\)に入り、(2)はu=1/2／2でもregularなので、
その全n2を含む中心も同じ閉和へ入る。
Pf(P_u^x)=Pf(P_u^y)=a0²/u≠0であり、translation frameのpoleをwallへ代入していない。

u=0／infinityはhighestにlinear Casimir、u=1は実quadratic pencilである。
新しい全linear-Hamiltonian／regular-pencil theoremにより、これらの全低次fiberも別に既知被覆される。
よってこの実M型は全period・全translation・全Cを含めて既知閉和へ入る。

## 3. 三つのmixed pole境界も全normal L零

旧全kernelのsquare条件から指定normal変数を零とした後、次の全式を得る。

| pole源 | 旧support帰結 | R | Pf nonzero constant | 残式のunit-forcing |
|---|---|---|---|---|
| x³y, full mixed | n0=n1=0 | 2a0+n2n3 | a0R/6 | n2R,n3R |
| x²q, q=xz+yw | n0=n1=0 | 2a0+n2n3 | −a0R/2 | n2R,n3R |
| x²y², full mixed | n2=0,n3=−n0 | a0+n0n1 | a0R | n0R,n1R |

各RはPf非零からunitである。従って前二行はn2=n3=0、最後はn0=n1=0。
保存されたnormal Lの全係数が零になり、全Cが残ったまま\(\mathcal B\)へ入る。
sourceの一部のquadratic係数だけからunitを推測せず、全非退化性を使っている。
残式の数値constant倍は非零であり、旧reportのQQ式を分母払して得る。

## 4. 全closed integrating-factor境界のrouting

両対数最高部閉包は旧proper quartic integrating-factor incidenceに入る。
部分分数の全極型を、既whole extension結果により次の通り戻す：

* reduced四独立line全support：六edge／五edge空、四cycleはregular-pencil、triangle/zeroはlinear-Casimir／零最高部。
* 実L112：q restriction rank二は空、rank一／零のaligned/transverseは空又はfullpullback。
  full mixed Mは§2。residue零／一致はlinear-Casimir／pencil、線形forms従属はfactor／pullback。
* x⁴はlinear-Casimir。
* x³y：rank二は空、aligned/transverseは空、Mは§3、全pullbackは別全定理。
* x²y²：rank二とaligned nonzero/transverse全parameterは空、Mは§3。
  aligned-zeroだけが§5の一軌道に残る。
* x²yz：独立全二源は空、全dependent配置はpullback。
* x²q：全rank restriction九源／唯一u零wallで空、fullmixedだけ§3。
  magnetic residue零／dependent primitive／pullback補集合は新全定理へ戻す。

factor最高部自体の全例外も `UNIMODULAR_ZERO_AND_FACTOR_TOP_ENTRANCE.md` と
全pullback membershipにより既知閉和へ戻せる。
上の空性は旧全低次sourceのPf恒等零を使い、generic topの空性を境界へ外挿しない。

## 5. 唯一残るx²y² aligned-zero源は一つのactual orbit

旧全残式では順にn0²=0、Pf非零からn2非零、a0*n2=0、2n2²n3=0となり、

```
n0=n3=a0=0, n1*n2≠0, Pf(P)=n1*n2³.
```

元の全normal tensorは正確に

```
Pxy=n2*x, Pxz=−x², Pxw=xz−n1*n2,
Pyz=xy−n2², Pyw=−yz+n2*(w−y),
Pzw=−xw+xy+n1*y+n2*z.                             (3)
```

全四translationを戻したsource以外に未処理のC/Lは無い。
最高部の全normalizerとtensor scaleを戻す前の二parameterに対し、actual diagonal map

```
(x',y',z',w')=(a*x, a*y/s², a*z/s, a*w/s²),
tensor scale=s, a*s≠0
```

は全六括弧を

```
n1'=a*n1, n2'=a*n2/s
```

へ送る。最高Q自体は全六成分で固定される。
a=1/n1、s=n2/n1を取ると全unit源がn1'=n2'=1へ移る。
このmapのdetはn1/n2⁵、inverseは

```
x=n1*x', y=n2²*y'/n1, z=n2*z', w=n2²*w'/n1,
inverse tensor scale=n1/n2.
```

Pfもscale²*detの共変則でexactに一へ移る。
従ってこの例外の全highest normalization、全translation、全lower源は、
次の**一つのactual PGL orbit**へ尽くされる：

```
P_*=(x,−x²,xz−1,xy−1,−yz+w−y,−xw+xy+y+z),
div P_*=0, Pf(P_*)=1.                              (4)
```

このone-orbit theoremは成分性又は既知親の同定ではない。
新しいambient点、接rank、minor、normal sliceをまだ一つも実行していない。
P_*の既知閉成分所属、又はそのorbit閉包を含む全成分の識別が、残る正確なE義務である。
単一canonically forced pointの判定以外に別pointを増やす理由は無い。
