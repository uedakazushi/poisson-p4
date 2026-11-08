# Stable line-conic：true T17 の quadric pencil を保つ actual Ga entrance

2026-11-08、PRIMARY E_multiple_hyperplane。純手の全 source 恒等式。
固定 H／固定 plane の normal-pair 診断を延長せず、moving H を保持する。
新 operator、point、kernel、minor、profile、jet、CAS は無い。
本稿は stable への actual endpoint 又は帰属を未証明として保つ。

## 1. True source の全 parameter

前稿 `NIL22_STABLE_LINE_CONIC_HQQ_MOVING_SOURCE_ENTRANCE.md` の completed T17 を

```
F=k*(AB+tau*H²), G=l*CD+nu*F,
PAB=F, PAC=sAC, PAD=rAD,
PBC=-sBC, PBD=-rBD, PCD=G, H row=0,
k=s+r
```

と書く。元 parameter は k=p1、tau=t0/k、nu=t1/k、l=p0+p2、
s=(p1-p2)/2、r=(p1+p2)/2。
実 source open は k 非零であり、tau,nu を独立の raw edge と取り違えない。

## 2. Quadric 二つを厳密に保つ Ga と全逆

任意 u に対し actual linear map

```
a=A, b=B-2*tau*u*H-tau*u²*A,
h=H+u*A, c=C, d=D                              (1)
```

を使う。det は一。inverse は同式の u を -u に変え、

```
A=a, H=h-u*a, B=b+2*tau*u*h-tau*u²*a,
C=c, D=d.
```

合成は Uu*Uv=U(u+v)。全二次式で
`AB+tau H² = ab+tau h²`、`CD=cd` なので F,G は個別に保存される。
これは polynomial coordinate 変換ではなく actual GL5。
任意の tau(epsilon),u(epsilon) に対しても、同一 parameter での正逆を持つ。

## 3. 全十 pushed bracket と physical H の符号

F=k*(ab+tau*h²)、G=l*cd+nu*F と読むと、(1) の全 push は正確に

```
{a,b}=F, {a,c}=sac, {a,d}=rad, {a,h}=0,
{b,c}=-s*c*(b+2*tau*u*h),
{b,d}=-r*d*(b+2*tau*u*h),
{b,h}=-u*F,
{c,d}=G, {c,h}=-u*s*a*c, {d,h}=-u*r*a*d.
```

省略又は独立 Taylor 化した source term はない。
例えば `{b,h}=u{B,A}=-uF`、
`{b,c}=-sBc-tau*u²*sAc=-sc(b+2tau*u*h)` から従う。

actual divisor は

```
Dnew=(h-u*a)*F*G.                               (2)
```

Poisson hyperplane は h=0 でなく **h-u*a=0**。
forward の h=H+uA と inverse の元 H=h-ua を混同しない。
(2) を (h+ua)FG に置き換えたり、新しい h 行を勝手に Casimir と呼ばない。
全 source の raw Euler 差は元 H-Casimir gauge の exact 差を同じ GL5 で移せば戻る。

## 4. 次の陽性義務

u が発散すると physical H は動くが、同時に full bracket の uF、usac、urad も動く。
したがって F,G と moving hyperplane の因子型一致だけでは regular tensor limit を得ない。
真 period k,s,r,l と edge k*tau,k*nu の coupled curve を先に指定し、
全括弧の pole cancellation と full stable six magnetic inverse を同時に示す必要がある。

stable line-conic 全源の二 magnetic invariant と全零境界は前稿 §4 のまま保持する。
この Ga の subgroup を何度も変えて新 frame/profile を走査することはしない。
具体的な endpoint formula が出た場合だけ ONE finite whole-variable readback の対象とする。
T16/T17 whole exclusion、T36/T37 不存在、新 component count は主張しない。
