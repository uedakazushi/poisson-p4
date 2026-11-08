# K22 の二 parameter 線：coisotropy と重み (0,0,1,2,2) の全 flag

2026-11-08、PRIMARY E、純手計算。CAS、profile、normal jet、追加点を用いない。
canonical K22 tensor 自体について証明する。K22 whole-source exhaustion/G上界を再証明せず、
K22 の既成分帰属も仮定しない。結論は以下の一つの line/flag weighted class だけであり、
全 K22 の plane 退化不存在を主張しない。

## 1. Actual tensor と線の全 normal bracket

座標 `(c,a,b,u,v)`、`P=Q+∂c∧W` とする。

```
Qab=-a*v+b*u, Qav=Qbu=uv, Qau=Qbv=Quv=0,
Wa=-2a²+4ab, Wb=4ab-2b²,
Wu=-2au+bu-av, Wv=-bu+av-2bv.
```

actual D は

```
F*=c*u²*v²+a²*v³+3abuv(u+v)+b²*u³.
```

`α,β` を任意とし `γ=-β²-3αβ`, `δ=-α²-3αβ` と置く。
線 L は `a=αu,b=βv,c=γu+δv`。F*|L=0 である。
normal 座標を

```
C=c-γu-δv, A=a-αu, B=b-βv
```

とする。元の全 tensor から直接展開すると

```
{A,B}=-v*A+u*B,
{C,A}=[-2αu+(α+4β)v]*A+3αu*B-2A²+4AB,
{C,B}=3βv*A+[(4α+β)u-2βv]*B+4AB-2B².                       (1)
```

従って三 normal bracket は L で全て零、すなわち L は coisotropic。
Euler representative を変えた線形式の bracket の差は ideal(L) に入るため、この性質は projective に不変。
線が D 内にあることだけから coisotropy を推論せず、(1) で確認した。

## 2. 一 normal を重み一、二 normal を重み二にする正しい flag 条件

projective normal 方向 `n=(C0,A0,B0)` を一つ選び、
`(C,A,B)=n*h+z*n1+w*n2` の定数可逆 normal frame を取る。
u,v の重み零、h の重み一、z,w の重み二、tensor scale=e² とする。

normal-pair {z,w} の L 上定数項は (1) により零。
その h 一次項だけは pullback で e^-1 の pole を持つので消す必要がある。
h² 項は e^0 であり、actual magnetic endpoint として**許容する**。
従って plane 自体の exact coisotropy を要求してはいけない。

三 normal 座標の bivector dual は `({A,B},-{C,B},{C,A})`。
その h 一次項と n の pairing を u,v ごとに零とすると、必要十分な全条件は

```
B0*[C0-(6α+β)A0+3αB0]=0,
A0*[-C0-3βA0+(α+6β)B0]=0.                                  (2)
```

`p=3α+2β`, `q=2α+3β` と書く。p*q 非零なら全 flag は正確に次の四つ。

```
(1,0,0), (-3α,0,1), (-3β,1,0),
(3α²+14αβ+3β², q, p).                                      (3)
```

実際 A0*B0 非零では両式から `p*A0=q*B0`、残る C0 は一意。
p=0,q非零なら第4はB0=0の既flagに合流し、q=0,p非零も対称。
p=q=0 はα=β=0で、全flagは `{C0=0}≅P1` と C軸の一点である。
これで全 zero/row-equality complement を保持する。

## 3. Actual quintic は全許容 flag で five hyperplanes

この PGL の determinant は e^5、tensor scale の二乗は e^4。
従って actual D は非零定数を除き `e^-1 F*(u,v,e*h,e²*z,e²*w)` である。
F*|L=0 だから、rank 四の端点は F* の h 一次項だけで決まり、

```
D_limit=h*u*v*q2(u,v),
q2=p*B0*u²+(C0+3βA0+3αB0)*u*v+q*A0*v².                    (4)
```

これは divisor のみの仮定ではなく、同じ actual tensor/PGL scale の Pfaffian 変換則である。
一般 n について normal→binary quadratic の写像が同型になる条件が p*q 非零である。
しかし pole-free の flag 条件 (2) を戻すと、(3) の四ケースで q2 はそれぞれ

```
u*v, p*u², q*v², (p*u+q*v)².                               (5)
```

に限られる。zero walls も (4) のまま戻せ、q2零なら rank 四 endpoint でない。
h² の magnetic 項を不当に零にしなくてもこの結論は残る。

従ってこの fixed finite-line/flag、一つの重み class の rank 四端点は、
五つの一次因子の積である。特に二超平面と既約 cubic、及び stable line-conic は得られない。
stable three-lines も得られない。前三ケースは三重以上の因子を持つ。
第4は h,u,v,pu+qv のうち u,v,pu+qv が同一点で三重に交わるが、
stable three-lines の四直線は三本ずつ独立である。

## 4. 射程と次の別幾何

上記は有限 α,β と固定 flag に対する通常の weighted contraction の全判定。
moving line/flag のより複雑な Laurent 依存、別の weights、Sing(D) 内の線、
他の K22/既知 parent の actual source 変形を否定しない。
K22 自身の known-MPS membership は別義務なので、得られた boundary 点が既知でも逆向きの E 被覆には使えない。

source 入力は `PLANE_CONCURRENT_M3_K22_ORE_NONDOMINANCE.md` §7 の明示 actual tensor だけ。
全源還元や G 上界を用いずに、この明示 tensor の線・flag・Pfaffian を手で監査した。
