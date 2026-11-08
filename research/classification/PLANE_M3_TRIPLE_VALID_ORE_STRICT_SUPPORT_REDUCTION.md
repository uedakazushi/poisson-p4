# generic triple の旧有効六 Ore 部分源：全 strict support の有限還元

2026-11-08、producer 純証明。全 generic triple の被覆とは区別する。
入力は `CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_ORE_20261108.json` の
**有効な六つの具体的 Q-family** と、その全 α wall の元 matrix。
旧自由-r coverage は sign erratum で撤回済みであり、本 note で復活させない。
正しい V2 被覆からこの有効部分源へ入る実写像は
`PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_ORE_ENTRANCE.md` §1 に別記した。
新 e二／η−一／r=dy/x 源はこの六つには含まない。

## 1. 全源と許される座標操作

`y=(a,b,z,w)`、`Π=Q(y)+∂c∧W(y)`。
各 Q は全 affine source の cohom 座標を持ち、W は weight-zero W0 に全十八 source を加える。
W の gauge は `Wz(z²)=Ww(w²)=0`。
元 fullM1 と FULLP2 九から Q,W への joint map は旧証明書に可逆行列として保存済み。
全 source を一方向の必要条件だけで減らしているのではない。

Jacobi の c 最高次数比較は strict `[Q,Q]=0` を要求する。
以下で F は normal binary quadratic、全 coefficient を最初に保持する。
base translations `a'=a+l(z,w), b'=b+m(z,w)` は plane を保つ actual linear change。
Qzw 零かつ Q の a,b 次数高々一なので、Qab の central F に対する作用は
その四つの translation coefficients に関して**正確に linear**である。
無限小 image を full orbit と取り違えない。inverse は l,m の negative。
この操作で W0 に加わるのは全 W十八に含まれる positive-weight 項なので、
W を部分的に固定することはない。
使用する diagonal changes も P0 の weight-zero source を保ち、
必要な非零係数を一にした後の W は依然として W0+全十八である。

共通 constant v が `v∧Q=L_vQ=0` を満たす全分岐は既 essential-three bound へ入る。
以下ではその分岐も除外せず明記する。

## 2. e 一／r=dx whole t

`d=3t+2` として全 Q は

```
Qab=−u d a z+F, Qaz=0, Qaw=−3u z²,
Qbz=u(6t+1)z², Qbw=9u t zw, Qzw=0.
W0=(a²,−2ab,3taz,(6t−1)aw−3bz).
```

全 strict Jacobi は恒等零。actual D3 は `3u(6t+1)abz³`。
従って u 非零、t≠−1/6 を要し、それ以外は元 actual leading image 零。
actual diagonal で u 一へ正規化する。
例えば `(a,b,c,z,w)` の scales `(A,B,C,Z,BZ/A)` と全 tensor scale A/C は
P0 を保ち、u を AB/(CZ) 倍にするので正規化は全非零 u を覆う。

base translation の central image は `span(z²,(3t−1)zw)`。
t≠1/3 は F=γw²、γ 零／非零の二 whole t-family。
t=1/3 は F=δzw+γw²。actual `w'=w+s z` は Qcohom と P0 を保ち、
γ 非零なら s で δ を消す。従って F=0,zw,w² の三型。
whole 二 family の特殊 fiber 0,w² に加え、F=zw の一 fixed を戻せば全 source を覆う。

## 3. e 二／r 零 whole t：全 strict support

`d=3t+2≠0` とし、cohom を A,U,V,R,S と書く。
旧全 Q-source に対する strict support ideal は正確に

```
d A U+3 R S=0,
(3t−1)U S−(3t+5)V R=0,
U V=0.
```

actual D3 は `A(3t−1)a²z³`。
t=1/3 の元全 matrix は actual D3 image 零なので、A 非零、t≠1/3 を残す。
U 非零なら V 零、次いで S 零、最後に dAU 非零となり矛盾。
従って U 零、RS 零、`(3t+5)VR=0`。
R 零の全 Q は b 非依存かつ全項に ∂b を持ち、essential-three 分岐。

R 非零なら S 零、t≠−5/3 では V 零。
actual diagonal で R=A=1 とすると、
`l=(6t−2)/(3t+2)` を用いて全 Q は

```
Qab=(l−1)bz+F, Qaz=l z², Qaw=(l+1)zw,
Qbz=0, Qbw=z², Qzw=0.
W0=(a²,−2ab,3taz,(6t+2)aw).
```

central image は z²,zw を含み、F=0,w² の二 whole t-family。
critical t=−5/3 で V 零なら同じ family の特殊 fiber。
V 非零なら actual diagonal で A=R=V=1 にでき、全 central F を translations で消す。
追加一 fixed は

```
(Qab,Qaz,Qaw,Qbz,Qbw)=(2aw+3bz,4z²,5zw,−4zw,z²−5w²),
W0=(a²,−2ab,−5az,−8aw).
```

### d 零の全 fiber

t=−2/3 は generic Q を代入せず、元全 matrix の八次元 Qsource を戻す。
その全 cohom を A,B,C,r,E とすると

```
Qab=−Abz−Bbw+F, Qaz=−Az²−Bzw, Qaw=−Azw−Bw²,
Qbz=−rz²−Czw, Qbw=−rzw−Cw²+Ez².
```

strict support は `Ar−BE=0, AC+rB=0, BC=0`、actual D3 は `−3Ea²z³`。
E 非零の下で B 非零なら C,r 零となり BE 零に矛盾、従って B 零。
A 零なら全 Q は共通 ∂b 分岐。
A 非零なら r=C 零、−A=E=1 へ正規化し

```
(Qab,Qaz,Qaw,Qbz,Qbw)=(bz+F,z²,zw,0,z²),
W0=(a²,−2ab,−2az,−2aw).
```

central image は z² のみ。w Borel が P0 と Qcohom を保つため F=0,zw,w² の全三型。
denominator d 零を単なる generic 例外として捨てない。

## 4. e 一 scalar／r=dlog x+dx fixed

元全 cohom u 非零を一へ正規化した Q は

```
(Qab,Qaz,Qaw,Qbz,Qbw)=(az+F,0,z²,−z²,z²−zw),
W0=(a²,−2ab,az,aw−3(a+b)z).
```

strict Jacobi は恒等零。central image は z² のみ。
w Borel は Qcohom/P0 を保つので F=0,zw,w² の三型で全源を覆う。

## 5. e 二 scalar／r=dlog x fixed

cohom A,B,C,D,E と F を用いて全 Q は

```
Qab=(2B+D)bz+2Dbw−Eaz+F,
Qaz=(2B+D)z²+2Dzw, Qaw=Bz²+(2B+2D)zw+2Dw²,
Qbz=−Az²+Ezw, Qbw=Cz²−Azw+Ew²+Ezw.
```

strict support は `DE=0, AD−BE=0, B(A+2E)−CD=0`。
leading amplitude は `L=A+C+E` 非零。
D 非零なら順に E,A,C 零となり L 零に矛盾。従って D 零。
B 零は共通 ∂b 分岐。B 非零なら E=A 零、C 非零。
B=C=1 へ正規化して

```
(Qab,Qaz,Qaw,Qbz,Qbw)=(2bz+F,2z²,z²+2zw,0,z²),
W0=(a²,−2ab,−2az,−2aw−3az).
```

central image は z² のみ、w Borel により F=0,zw,w² の全三型。

## 6. R3L／φ零／r=dx fixed

cohom A,B,C,D,E,T と F を用い、全 strict support は

```
E(D+E)=0,
E(A+T)−B(2D+3E)=0,
B(A+2T−2B)−CE=0.
```

leading amplitude `L=A−B+T` 非零。
E 非零なら D=−E、次の式は EL 零となるので矛盾。
従って E 零、BD 零、`B(A+2T−2B)=0`。
B 零は共通 ∂b 分岐。
B 非零なら D 零、A=2B−2T、L=B−T 非零。
actual `b'=b+s a` は P0 を保ち C を C+sB へ変えるので C 零にできる。
L 一へ正規化し `ρ=B/L∈C*` を全て残すと

```
(Qab,Qaz,Qaw,Qbz,Qbw)=((ρ+1)az+F,0,ρz²,−2z²,(ρ−3)zw),
W0=(0,−3a²,0,−3az).
```

central image は z² と `(ρ−1)zw`。
ρ≠1 では F=0,w² の二 whole ρ-family。
ρ=1 では w Borel が P0 と Qcohom を保つので F=0,zw,w²、追加一 fixed F=zw を戻す。
ρ零は前の B零分岐、L零は actual D3零であり、除いたままではない。

## 7. R3L／φ=x／r零 fixed

旧全 Qsource の strict reduced support は
`Q1²=0, Q0 Q1−Q3²=0, 2Q1Q3+Q1Q4−Q2Q3=0`。
従って全 reduced support で Q1=Q3=0。
残全 Q は b 非依存、全項に ∂b を持つので essential-three 分岐。
nilpotence exponent 又は非被約 scheme の分類は不要である。

## 8. 次の元方程式判定の有限対象

共通 constant-v 分岐を除いた有限対象は

* §2 二 whole＋一 fixed。
* §3 二 whole＋critical V 一 fixed＋d零三 fixed。
* §4 三 fixed。
* §5 三 fixed。
* §6 二 whole＋一 fixed。

計 **六 whole parameter family＋十二 fixed＝十八**。
全 W十八は各 source に必ず戻す。
その full original Jacobi affine W system と actual D の proper factor 又は source上界が
証明されるまでは、この十八 family の全 G閉鎖は未証明。
さらに本 note 外の新 V2 fixed source も残るため、全 generic triple の完了は主張しない。
