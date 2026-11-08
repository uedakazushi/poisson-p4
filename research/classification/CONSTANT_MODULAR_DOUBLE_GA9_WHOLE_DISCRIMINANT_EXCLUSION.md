# constant DOUBLE は全 Ga9 cuspidal quartic 成分の外：D5 の全変形による排除

2026-11-08、PRIMARY E。担当 `/root/plane_collision_review`。
**純全成分証明の producer。独立査読待ち。** 新 CAS、point/operator/minor、profile 又は jet 延長無し。
Ga9の一標本の D100 を全sourceへ外挿しない。
actual quintic の因子と plane quartic の全局所変形だけを使う。

## 1. 対象と正しい Ga9 全成分入力

受理済み DOUBLE canonical の actual degeneracy divisor は

```
D0=h²*C, C=x²y+xh²+h³.
```

Cは絶対既約で、D0は(h,x,y)の三つを本質的に使う。
Ga9 cuspidal quartic **24成分**の稠密sourceは、全六 magnetic、全PGL frameと tensor scalarを保って

```
Da=h*Fa,
Fa=(h*w+6*y²)²−144*y³*(h+a*y).
```

`CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json` の
`actual_family_power3.plane_quartic` そのものを使う。
power1全源との可逆全150移送は既受理で、power3の稠密源/PGL閉包が全成分である。
従って本稿でaの例外値を別未知源としない。
磁気はDを変えず、全六係数を任意に残す。

a≠0、1−4a≠0の稠密開でFaは既約で、二つの相異なる singularity を持つ。
q=[h:y:w]=[1:0:0] のh=1 chartではW=w+6y²として

```
Fa=W²−144*y³*(1+a*y),
```

よって q は A2 cusp。
p=[0:0:1] のw=1 chartでは H=h+6y²−72y³ と置くと

```
Fa=H²−144*a*y⁴+864*y⁵−5184*y⁶,
```

よって a≠0 の全域で p は A3 tacnode。
A3,A2 の存在だけを使い、他のsingularitiesが無いことは本証明に不要。

## 2. 任意の moving-frame 境界から actual P2 の変形へ

仮にDOUBLEがwhole Ga9閉成分に属すとする。
稠密constructible源からの curve selection と有限base changeにより、
非零parameterで§1の稠密sourceの actual PGL像、零でD0となる一parameter曲線を取れる。
Pfaffian mapは非零D0の近傍でregularなので、actual quintic の曲線も得る。

各非零fiberのessential linear subspace Wt⊂V*はdimension三。
Gr(3,V*)のpropernessで Wt はrank三bundle W0 へ延長し、D0∈Sym5(W0)。
D0のessential rankも三だから W0はちょうどspan(h,x,y)。
DVR上のGrass bundleの局所frameを選べば、regularかつ中心可逆なGL5変更により
全Wtを一つの三次元Wへ戻せる。
これは任意のsingular PGL変換を中心でも可逆と仮定する操作ではない。
先にそのimage Wtだけをproperに延長し、そのbundleのregularframeを選んでいる。

因子空間P(W)×P(Sym4 W)のproper multiplicationにより、
markされたlinear factor htとquartic Ftも延長する。
D0の唯一のlinear因子はhなので h0=h、従って

```
F0=h*C.
```

因子htがWt内にあることは、Wtをannihilateする定数微分dが
0=d(htFt)=d(ht)Ft+ht*dFt を満たし、ht∤Ftより d(ht)=0となることからも分かる。
これで全moving-frame境界を実際の plane quartic Ft→hC に還元した。
magnetic/unit/periodの消失を勝手に固定していない。

## 3. 中心 quartic の全 singularity

F0=h(x²y+xh²+h³) の singularity は H∩C の二点だけである。
Cの唯一のsingular pointは p0=[h:x:y]=[0:0:1]。
H∩C は p0 と q0=[0:1:0] である。

q0ではx=1として F0=h(y+h²+h³)、Hessianは非退化なので A1 node。
p0ではy=1、X=x+h²/2 と置くと

```
F0=h*(X²+h³*(1−h/4)).
```

z=h*(1−h/4)^(1/3)、w=X という可逆な収束座標とunit倍により

```
z*(w²+z³)
```

となる。従ってp0は厳密にD5である。
ここでの非線形解析座標は局所singularityの同定だけに使い、actual PGLと呼んでいない。

§1のA3点とA2点はP2のpropernessで中心へ延長する。
中心のsmooth点の近傍にsingular pointは現れない。
q0のHessian非退化性はその近傍の全singular pointsにも開条件として保たれるので、
A3又はA2はq0へ極限出来ない。
従って二つの相異なるA3,A2 marksはともにp0へ集まらなければならない。

## 4. 使用する primary miniversal theorem と全 T1 基底

一次資料：
[Fedorchuk, arXiv:1007.4828v2](https://arxiv.org/pdf/1007.4828v2) §1.2.2 printed7、
§6 Proposition6.1 printed29–30。
D5の全miniversal familyは、符号をcomplex unitで調整して

```
F(z,w)=z*w²+b*w+R4(z),
R4(z)=z⁴+a3*z³+a2*z²+a1*z+a0.                    (1)
```

これは一つのchosen sliceでなく、**全局所変形**をpullbackするfamily。
収束analytic範疇での isolated complete-intersection semiuniversal theoremは
[Tyurina, Math. USSR-Izv.3(1969), Theorem8.1, printed997](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=2190&what=fullteng)
を直接確認した。CIの場合にnormalityを仮定しない。

本例で基底も手で確認出来る。f=z*w²+z⁴のJacobian/Tjurina商は
zw=0、w²=−4z³、z⁴=0を満たし、基底は

```
1,z,z²,z³,w.
```

従って(1)の五parameterは全T1方向を持ち、Fedorchukの明示的versalityと一致する。
この外部定理以外に、A/D adjacency table又は未読の分類結果は使わない。
次節は(1)の**全parameter**についての新しい初等的検算である。

## 5. D5のどの小変形にも A3+A2 は同時に存在しない

### z非零領域

z倍して平方完成すれば

```
z*F=(z*w+b/2)²+P5(z), P5=z*R4(z)−b²/4.
```

z≠0では(z,w)→(z,z*w+b/2)は可逆で、unit z倍もsingularityを変えない。
P5のmultiplicity mの根は、局所equation W²+(z−z0)^m*unit に対応する。
A3にはm=4、A2にはm=3が必要。
相異なる二点は異なるz0を持つ（critical W=0でwは一意）ので、
両者がz≠0にあれば degree五 polynomial に合計七重根が必要となり不可能。

### z零領域：全補集合

z=0のsingular pointがあれば Fw=b よりb=0、F=0よりa0=0。
その時

```
F=z*(w²+R3(z)), R3=z³+a3*z²+a2*z+a1.
```

* a1≠0：z=0上の二点 w²=−a1 は、二因子が横断するA1だけ。
  A2/A3はここになく、z非零領域のdegree五の議論が適用される。
* a1=0、a2≠0：z=0の唯一点は、smooth因子 z と smooth因子
  w²+a2z+… のintersection multiplicity二、従ってA3。
  他方z≠0のA2にはR3の非零三重根が必要だが、R3はdegree三で既に零根を持つ。
  よって非零三重根は無い。A3+A2は不可能。
* a1=a2=0：唯一のz=0点はmultiplicity三なので、multiplicity二のA2/A3ではない。
  R3=z²(z+a3)の非零根はあってもsimpleであり、z非零にもA2/A3はない。

b又はa0非零、a1非零、a1零a2非零、a1a2零の全補集合を尽くした。
a3の値、全zero fiber、reducible local fibersも捨てていない。
従って(1)の全fiberに相異なる A3 と A2 が同時に現れることは無い。
semiuniversal theoremにより、同じ結論は任意のD5局所変形に成立する。

## 6. 全 Ga9 成分排除と残余

§3はA3+A2の両marksがD5の一近傍へ集まることを強制するが、§5に反する。
従ってD0は、全Ga9 componentのactual divisor image closureにも入らない。
よってDOUBLE canonicalのactual orbitは **whole Ga9 cuspidal quartic24成分に入らない**。

全六magnetic・全PGL・aのzero/例外/infinity境界を含む閉成分の排除であり、
a1 witnessのambient rankを一般点や全成分へ外挿していない。
従来のwholeT36/T38排除とknown24 filtersに加えると、known24親はT37だけが残る。
known23の全監査、unknown上界、新成分性又はT37 membershipは別義務のままである。
