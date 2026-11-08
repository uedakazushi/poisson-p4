# Stable line-conic：T16/T17 の真 HQQ 源と moving-source に残る tensor 条件

2026-11-08。PRIMARY E_multiple_hyperplane、純手の exact entrance。
新 point/operator/kernel/minor/profile/jet/CAS は実行しない。
**本稿は陽性帰属も whole-parent 陰性も証明しない。**
既知 T16/T17 はどちらも projective dimension24 なので、全 stable source dimension23 だけでは除外できない。
ID の数値範囲から component dimension を推測しない。

## 1. 元全 source と全 H-Casimir gauge

source 順を (A,B,C,D,H)=(x0,x1,x2,x3,x4) とする。
T16 は公開 raw catalogue source、T17 は correction t0*t1*H²/p1 を含む true completed source を使う。
両者とも X_i=Praw_iH/H とし P_H=Praw+E wedge X と置く。
H 行を消すこの全十括弧の exact Euler 変更は、元代表へ同じ項を引けば戻る。

T16 は

```
PAB=F=kAB+tH², PCD=G=lCD+uH²,
PAC=-sAC, PAD=sAD, PBC=sBC, PBD=-sBD;
k=(p2-p1)/2, s=(p1+p2)/4, l=p0+2s, t=t0, u=t1.
p1=2s-k, p2=2s+k, p0=l-2s.
```

cross の Pfaffian 寄与は相殺し、actual D=HFG。
raw Euler X=(0,kB,-sC,sD,0)。

T17 は

```
PAB=F=kAB+tH², PCD=G=lCD+(u/k)F,
PAC=sAC, PAD=rAD, PBC=-sBC, PBD=-rBD;
k=p1, s=(p1-p2)/2, r=(p1+p2)/2, l=p0+p2, t=t0, u=t1.
k=s+r, p2=r-s, p0=l-r+s.
```

実 completion の denominator k は generic source open で非零。
raw Euler X=(0,kB,sC,rD,0)。actual D=HFG。
特に u と (u/k)F の全 parameter compatibility を別の自由係数に緩めない。

## 2. Stable の H=double-h assignment に対する exact divisor match

target base は (h,x,y) で、全六磁気に共通して

```
D=const*h²*L*Q,
L=x-y-12h, Q=36h²+3h(x-y)+2xy.
```

H は target h として次の literal 因子対応がある。

T16:

```
A=h, B=L, C=x-3h/2, D=y+3h/2,
B=C-D-9H,
F=hL, G=2CD+(81/2)H²=Q.
```

T17:

```
A=h, B=L, C=x-6h, D=y+6h,
B=C-D,
F=hL, G=2CD-9F=Q.
```

したがって divisor の次数・rank・marked contact は両 true source へ exact に対応し得る。
これだけで Poisson tensor を一致させたとは呼ばない。
source の三二次因子と source period の可逆性は §1 のまま保つ。

## 3. 同じ固定 plane での normal-pair 誤差

上の target base plane へ source を制限したとき、normal を次の二つに取る。

T16: f=A-H, g=B-C+D+9H。
source 全 bracket から厳密に

```
{f,g}|S = H*((k+s)C+(s-k)D+(t-9k)H).
```

T17: f=A-H, g=B-C+D。
全 bracket から

```
{f,g}|S = H*(rC-sD+tH).
```

どちらも generic k 非零では identically zero ではない。
固定 plane に二 normal を一様に epsilon 倍し tensor を epsilon 倍するだけでは
normal-pair に epsilon^-1 の実 pole が残る。
この指定 class の leading coisotropy を課すと T16 は k=s=t=0、T17 は r=s=t=0、従って k=0。
これは一つの anchor も失わせるため、factor matching に同じ generic source を凍結して使えない。

この診断は true moving-period family、別の degeneration flag、meromorphic base frame を排除しない。
特に k=0 source boundary の (u/k)F は個別の項を落とさず、全 tensor の正則性を先に証明する必要がある。
固定 class の次数・profile を増やさない。

## 4. 全磁気と未解決義務

stable full source は任意 R=m0h²+m1hx+m2hy+m3x²+m4xy+m5y² を持つ。
既受理全 polynomial fiber-translation normalizer は rank4 で、残る二つの不変量は

```
I4=m4-2m1/3-8m3,
I5=m5-m1/12-m2/12-m3.
```

全 normalizer と逆は parameter を反転しない。I4/I5 の零境界も元全六源も対象である。
anchor Borel action だけ、divisor の一致だけ、同じ二不変量の一点だけから whole E を推測しない。
残る named T16/T17 route の最弱陽性目標は、真 period/edge compatibility を保つ actual curve と
全二不変量（あるいは全六磁気）の inverse を同時に示すこと。
T36/T37 及び他未排除の global parent は本稿の陰性対象ではない。

T03 は別の actual curve により stable_three_lines の全六磁気源を被覆した。
この結果を stable_line_conic へ延長しない。
