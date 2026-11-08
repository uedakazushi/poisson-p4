# generic triple・非被約Rの全first-symbol有限正常形

2026-11-08、producer純幾何。全first connectionの分類であり、globalK/全higher/G閉鎖は未了。
前 `PLANE_NONREDUCED_M3_TRIPLE_ROOT_FIRST_CONNECTIONS.md` の入口を用いる。
以下の全offdiagonal sourceはnormal Borelとactual affine base stabilizerで尽くす。

## 1. R2L全八dim covariant source

f=a²b、h=a^e b^(2−e)、e=0,1,2、G3=h z³。
x=b/a,y=c/a、η=4/3−e−αと置く。

`p=(dlogf)/3+αdlog(a/b)+γdx+δdy, q=2p−dlog(f/h)`。

q−p=ηdlogx+γdx+δdy。
全global rは

```
r=(Pdx+Qdy)/x,
P=A+Bx+Cy+Dxy+Ey²,
Q=F+Gx+Hy−Dx²−Exy.
```

H0Ω1(3)の八係数を全て戻したもので、scale3はrの任意係数へ吸収している。
flatnessは一つの全polynomial identity

`x(Q_x−P_y)+(η−1)Q+γxQ−δxP=0`。

normal Borelはr→r−s(ηdlogx+γdx+δdy)、normal diagonal ratioはr全体の非零scale。
base stabilizerはx→λx,y→sy+tx+uで、使うγ/δ正常形を保つ部分だけを用いる。

δ非零源は `PLANE_NONREDUCED_M3_TRIPLE_IRREGULAR_C_EXCLUSION.md` により全空である。
以下ではγ一/δ零又はγδ零だけを残す。

## 2. γ一：全二whole lineとη一の一追加型

全coefficient比較からη≠1では

`r=A dlogx+B dx`。

normal BorelでBを消し、残A−ηBはnormal ratioで零/一へ。
従ってr=0,dlogxの二whole α lineが全ηで有効。
η=1では全sourceは

`r=A dlogx+Bdx+F(y dx+dy)/x`。

normal BorelでBを消し、F非零ならy translationで残dlogx係数を消す。
従ってr=(y dx+dy)/xを一fixed追加すれば全resonanceも尽くされる。
このη=1は各eについてα=1/3−e。
全eで六whole line＋三fixed。

## 3. log-only generic η

γδ零ではcoefficient identityは正確に

```
C=ηG,
(η−1)F=(η−1)H=0,
(η+2)D=(η+2)E=0.
```

η≠1,−2ではr=A dlogx+Bdx+G(ηy dlogx+dy)。
η≠0ならnormal BorelでAを消す。
G非零でη≠−1ならy→y+txでBを消せる。
従ってη∉{−2,−1,0,1}では

`r=0, dx, ηy dlogx+dy`

の三whole α sourceが全てを覆う。

## 4. log-only四resonanceの全source

η=−1ではnormal BorelでAを消せるがy shearがBを動かさない。
全rは

`0, dx, −y dlogx+dy, −y dlogx+dy+dx`

の四fixed。normal ratio、x/y非零scalesで残B/Gの全零/非零を戻す。

η=0ではBorel imageが零なのでAを落とさない。
全r=A dlogx+Bdx+Gdyに対して、G非零ではy shearでBを消す。
全六fixedは

`0, dx, dlogx, dlogx+dx, dy, dlogx+dy`。

η=1ではd(xr)=0。normal Borelによりprimitive ψをmodulo C*x及びconstantで考えられる。
全r=dψ/xについてψの八型を取る：

`y², y²+x², xy, xy+y, 0, x², y, x²+y`。

証明：y²係数非零ではy translation/shearでxy,yを消し、x²残係数零/非零の二型。
y²零/xy非零ではx²を消し、linear y零/非零の二型。
両quadratic係数零ではx²とyの各零/非零の四型。
全transformはx→λx,y→sy+tx+uで、normal ratioも自由である。

η=−2ではd(r/x²)=0。normal Borelはprimitive numeratorのconstantを変える。
全r=x²d(ψ/x²)についてψをmodulo C*x²及びconstantで考え、七型を取る：

`y², y²+x, xy, xy+y, y, 0, x`。

y²係数非零ならy translation/shearでy/xyを消し、残x零/非零の二型。
y²零/xy非零ならxを消し、残y零/非零の二型。
両quadratic係数零ならy非零の一型、y零/x非零又は零の二型。
同じactual groupとnormal ratioにより全inverseがある。

全resonanceはη値によりα=4/3−e−ηへ戻す。
全eでlog-onlyは九whole line＋七十五fixed。
各whole lineから四η値を除くが、その全fibersを上記二十五fixedへ復活している。

## 5. R3Lの残十fixedと全数

f=a³,h=a²では前noteの全potential φ五型を使う。
φ=x²+y,xyはq(C)の定数vector contractionがinjectiveなので全空。
残φ零の全r=dψ五型、φ=xのr0/xdx/ydx+dy三型、φ=x²のr0/dx二型、全十fixedを保持する。
p=dloga+dφ,q=dloga+2dφ、全r normalizerは前noteに明記。

よって全generic triple first connectionsは

- R2L γ一：六whole line＋三fixed。
- R2L log-only：九whole line＋七十五fixed。
- R3L：十fixed。

合計十五whole parameter source＋八十八fixedで覆われる。
有限source列挙は全r/sourceを含む構造分類であり、標本値のprobeではない。

## 6. source dimensionと次の弱い計算

R2Lのfixedplane全sourceは、orderedR四＋πS amplitude一＋constant normal root一の六に、
p-source dimensionと元r-solution stratum dimensionを加えた上界で数えられる。
γ一genericはp二＋r二でB≤10、η一でもp一＋r三で同じB≤10。
log-only genericはp一＋r三でB≤10、η−1/0 fixedはB≤9、η1/−2はr源五を戻してB≤11。
canonical Borel gaugeで消したr座標をsource dimensionから無断で除かない。
R3Lの十型は前noteのBを使う。

次の有限計算ではrawglobalK/P0/full30/9/actualD3lineに加え、
fullP1kernelのbad normal三行が零なら直ちにz因子で空とする。
fullD3像零又はB+k1+k2+5≤24も弱い十分条件である。
未達なら全sourceを保存し、元Jacobi/Ore/actualDを使う。
本finite化だけでgeneric triple又は全m三Gを閉じたとはしない。
