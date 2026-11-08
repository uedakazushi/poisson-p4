# constant DOUBLEのcubic成分：一般rank二leafは代数面に入らない

2026-11-08、PRIMARY E_multiple_hyperplane、別作者 `/root/quartic_sources` の純手証明。新CAS、point、rank/operator、局所jet又は算術実行は無い。rootの印刷式を独立に再導出した。**本稿はtarget上のno-leafを証明する。whole T37に必要なproper leaf-incidence定理は別作者作成中の依存であり、まだ本稿だけでT37排除又は新成分を主張しない。**

## 1. 全targetとcubic成分上の正確な括弧

全最高double-line sourceの独立有限還元で強制されたcanonicalを使う。h=1、座標(t,x,y,r)で

```
Ptx=1+x-2xy, Pty=4y²-y, Ptr=(1-6y)r+x,
Pxy=0, Pxr=-x², Pyr=1+2xy.
```

actualD=h²*C、C=x²y+xh²+h³。Cはyについてprimitive一次式であり、x²とh²(x+h)が互いに素だから複素数体上でも既約である。h=1上のCはx=0を含まず、全域で

```
x unit, y=-(x+1)/x².
```

従ってこのcubicのh非零開はSpec C[x,x^-1,t,r]と同型で滑らかである。P sharp dC=C partial_tを直接展開出来るので、CのidealはPoisson invariant。制限括弧は

```
{t,x}=A=(x+1)(x+2)/x,
{x,r}=-x²,
{t,r}=C0*r+x,
C0=(x²+6x+6)/x².
```

xがunitなので、rankはこの開全体で正確に二である。

## 2. 全kernel一形式と代数的座標変更

volume dt wedge dx wedge drに対するkernelは

```
alpha=-x²*dt-(C0*r+x)*dx+A*dr.
```

mu=A/x²=(x+1)(x+2)/x³と置くと

```
mu'=-(x²+6x+6)/x⁴=-C0/x².
```

u=t-mu*rはこの開で正則な代数的座標で、逆t=u+mu*rも同じring上にある。よって

```
alpha/(-x²)=du+dx/x.                              (1)
```

従ってrank二distributionは全域で

```
partial_r, D=x*partial_x-partial_u
```

により張られる。logの解析的枝をactual PGL又は代数的first integralと取り違えていない。

## 3. Invariant代数面の不存在

この滑らかな三次元開のcoordinate ringはUFD C[x,x^-1,u,r]。distributionにinvariantなirreducible代数面があれば、そのprime idealは非unit irreducible f一つで生成される。

partial_r fが(f)に入ることとr次数を比較すると、partial_r f=0。したがってfはC[x,x^-1,u]に属する。次にDf=g*fと書く。Df=0なら、xごとのpolynomial係数を比較してfは定数となるので非unitに反する。

Df≠0の場合、Dはxの最大／最小指数を広げず、u次数も増やさない。積の最大／最小指数はdomain C[u]上で加法的だから、gの最大x指数は高々零、最小指数は少なくとも零となる。gはxに依らない。u次数も比較してgは定数cである。

`f=sum_j x^j*p_j(u)` と書けば

```
p'_j=(j-c)*p_j.
```

非零polynomial p_jが存在するにはj=c、p_j定数が必要であり、他のjは全て零。従ってfはnonzero scalar*x^j、すなわちLaurent ringのunitとなる。矛盾。

よって**このcubic成分のh非零開には、rank二distributionにinvariantな代数面は一つも存在しない**。一般symplectic leafはそのcubic三次元成分にZariski denseである。特に、cubicの一般点を含みそのrank二leafを含むprojective代数面の閉包も存在しない。

## 4. Proper leaf flagに必要な正確な射程

T37側で単に「低次数の代数leafが一つ存在する」と示すだけでは、本targetを除けない。special actualDはdouble HとCの和なので、その一つのsurfaceだけならH内へ退化して逃げ得る。必要なのは、generic irreducible quinticの一般点を被覆するalgebraic leaf-surfaceをmarkし、targetのC一般点を含むflat limitを得ることである。

別作者が現在検討している方針は、非零quinticのDVR flat family、target C一般点を通るpoint lift、generic点を含むbounded-degree leafの選択、proper Hilbert scheme上のflat ideal closureとHamiltonian derivation invarianceの延長である。本稿ではそのT37全source側の構成、degree bound又は一般点被覆を先取りしない。

この方針が証明された後には、flat limitにembedded成分が現れても、§3との矛盾を避けられない。char零のderivationにstableなidealのradicalとその最小primeはstableである。したがってrank二のregular点を含むreduced supportのirreducible componentは次元二以上を持つ。flat surfaceのsupportにその点が残り、全surfaceがspecial D内にあるなら、h非零点を含む二次元componentはC内のinvariant代数面になる。零次元embedded点だけをそのrank二点へ残すことでは済まない。

以上のtarget側計算は別作者として受理する。T37全proper-incidence依存、global componenthood、prior比較と新成分数は依然別義務である。同一forced orbitのquadratic余核零結果も変更せず、その後のcubic／高次local計算を再開しない。
