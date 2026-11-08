# Stable line-conic 全六磁気源：true T17 の正しい quadric assignment と全 period 曲線

2026-11-08、PRIMARY E_multiple_hyperplane。新しい pure actual 全源曲線。
新 point/operator/kernel/minor/profile/jet/CAS は用いない。
以下は全変数恒等式を示す純証明であり、同一式の ONE finite readback と独立査読を後続義務とする。
全 stable_line_conic 六磁気源、PGL 飽和とその閉包が named T17 に入ることが目標。
新 component count は零。

## 1. 真 source と従来診断の assignment の違い

`reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T17.json` の全 source を使う。
座標順 (A,B,C,D,H)=(x0,x1,x2,x3,x4)。H-Casimir gauge の全六式は

```
F=kAB+tH², G=lCD+nu*F,
PAB=F, PAC=sAC, PAD=rAD,
PBC=-sBC, PBD=-rBD, PCD=G,
k=s+r.
```

raw period 正逆は

```
p1=k, p2=r-s, p0=l-r+s, t0=t, t1=nu*k;
k=p1, s=(p1-p2)/2, r=(p1+p2)/2, l=p0+p2.
```

raw Euler vector は X=(0,kB,sC,rD,0)、P_H=Praw+E wedge X。
完成項 t0*t1*H²/p1 は t*nu*H² と厳密に相殺する。
actual divisor は D=HFG。

前の固定 plane 診断は F→hL,G→Q という assignment のみであった。
今回は **F→Q、G→hL** を使うのでその陰性を反復しない。

## 2. Exact coisotropic plane と元 stable anchor

a=9/2 とし、source 中心を

```
s=r=l=1, k=2, t=2a², nu=-1/2
```

に置く。plane S を C=A-aH,D=B+aH とする。
f=C-A+aH、g=D-B-aH なら

```
{f,g}=G+F-sBC-rAD,
{f,g}|S=0,
F|S=2(AB+a²H²), G|S=aH(A-B-2aH).
```

normal だけを epsilon 倍し tensor を epsilon 倍する。
U=f/epsilon,V=g/epsilon、base は(H,A,B)。pole は E0={f,g}|S=0 で消える。
端点全 bracket は

```
{U,A}=-A(A-aH), {U,B}=-AB-aBH-2a²H²,
{V,A}=AB-aAH+2a²H², {V,B}=B(B+aH),
{U,V}=aH(U-V), {A,B}=0, H row=0.
```

actual constant frame

```
h=H, x=A+3H/2, y=B-3H/2,
Q1=U/3, Q2=-V/3
```

で四 anchor は元 stable_line_conic に全一致する。
q-pair は -3h(Q1+Q2)/2。
F|S=36h²+3h(x-y)+2xy=Q、A-B-2aH=x-y-12h=L。
従って D は元 h²LQ 型で、double h は Q を二点で切り、simple L は Q に接する。

## 3. 真 period derivative による四つの磁気係数

自由 S,R,T,V0 に対し同じ epsilon 上で

```
s=1+epsilon*S, r=1+epsilon*R, l=1,
k=s+r, t=2a²+epsilon*T, nu=-1/2+epsilon*V0.
```

normal plane 上の誤差は正確に

```
E=(l+nu*k)AB+a(l-r)AH+a(s-l)BH+((1+nu)t-la²)H²
 =epsilon*E1+epsilon²*E2,
E1=(2V0-(S+R)/2)AB-aR*AH+aS*BH+(T/2+2a²V0)H²,
E2=V0*(S+R)AB+V0*T*H².
```

したがって full normal-pair は E/epsilon を含むだけで全負冪が消える。
他の全括弧も normal の一次 inverse しか持たず tensor epsilon で相殺する。
これは true source 全 parameter の実曲線であり、period と edge の積を自由 Taylor 係数へ緩めない。
端点の追加磁気は -E1/9。四係数 map は可逆。

## 4. 任意全六磁気の explicit binding と inverse

元 Rmag=m0h²+m1hx+m2hy+m3x²+m4xy+m5y² に対し

```
c0=m0-9m3-9m5, c1=m1+12m3, c2=m2-12m5,
rAB=m4,
rAH=c1-3m4/2, rBH=c2+3m4/2,
rHH=c0+3c1/2-3c2/2-9m4/4;
S=-2*rBH, R=2*rAH,
V0=-9*rAB/2+(rAH-rBH)/2,
T=-18*rHH-81*V0.
```

すると -E1/9 は c0h²+c1hx+c2hy+m4xy になる。
最後に actual finite shift

```
Q1final=Q1+3m5*y, Q2final=Q2-3m3*x
```

を行う。元 q-pair の linear-Q 項も同時に置換すると磁気変化は

```
m3*x²+m5*y²-12m3*hx+12m5*hy+9(m3+m5)h².
```

よって全六磁気へ厳密に戻る。
全 map m→(S,R,T,V0,3m5,-3m3) の determinant は -2916。
逆は shift を (j1,j2) と書いて

```
m5=j1/3, m3=-j2/3,
rAB=-2V0/9+(S+R)/18, rAH=R/2, rBH=-S/2,
rHH=-T/18-9V0/2,
m4=rAB,
m1=rAH-12m3+3m4/2,
m2=rBH+12m5-3m4/2,
m0=rHH-3m1/2+3m2/2-9m3-9m5+9m4/4.
```

磁気係数又は invariant を一つも反転せず、全零境界を含む。

## 5. 全 actual GL5 と closed-component 論理

source→final target は

```
h=H, x=A+3H/2, y=B-3H/2,
Q1=(C-A+aH)/(3epsilon)+3m5*y,
Q2=-(D-B-aH)/(3epsilon)-3m3*x.
```

inverse は

```
A=x-3h/2, B=y+3h/2, H=h,
C=x-6h+3epsilon*(Q1-3m5*y),
D=y+6h-3epsilon*(Q2+3m3*x).
```

det=-1/(9epsilon²)、tensor scale=epsilon。
全 signed Pfaffian の変換係数は epsilon² det=-1/9。
端点 actual D は -h²*L*Q/2 で元全源と一致する。

任意に固定した六 m に対し、epsilon 非零で k,t,nu の非零 open は中心で
k=2,t=81/2,nu=-1/2 なので非空。
source biresidue Pfaffian open は p1(p0+p2)=k*l=k で同じく中心非零。
従って generic epsilon は true T17 source の actual PGL 像と非零 tensor 倍。
全括弧が epsilon 零まで多項式かつ端点 D 非零なので、T17 の閉性から全六磁気源、
全 PGL 飽和、その閉包が T17 に入る。
境界帰属の逆推論、local uniqueness、smoothness、ambient component count を用いない。

有限 readback は raw Euler 全150正逆、元 full6magnetic endpoint、
全strict10／projective350、全126 Pf、全 frame／parameter inverse を同時確認する。
この純 note はその producer 実行と別作者 review に先行して凍結する。
