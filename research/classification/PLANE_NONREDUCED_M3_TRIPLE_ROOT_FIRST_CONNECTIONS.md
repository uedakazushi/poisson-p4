# 非被約R・m三generic triple root：全first connectionとR3L十二型

2026-11-07、producer純幾何。rank二plane、actual reduced irreducible quintic、位数三。
generic normal leading cubicがtriple rootならUFDからG3=h z³、h quadratic、
primitive rootのbase degreeは3d≤2により零。div(h)の支持はR。
本noteはglobalK/higher-source/G閉鎖をまだ主張しない。

## 1. 全接続

zを保存する全conormal β=[[p,0],[r,q]]に対し

`2p−q=dlog(f/h), dp=0, dr+(q−p)∧r=0`。

pの全nonlog部分を残す。rはglobal f*r polynomialというH0Ω1(3)八次元源に属する。
triple rootだからと三distinct-rootの接続一意性を使わない。
normal Borel w→w+s zはr→r+s(p−q)、normal diagonal ratioはrの非零scaleを変える。
全inverseはactual constant normal GL2である。

## 2. R=3L：全potentialと全offdiagonal source

f=a³、h=a²。p=dloga+dφ、q=dloga+2dφ。
a非零chartのx=b/a,y=c/aにおいてφはdegree≤2 polynomial modulo constant。
actual affine base GL2とtranslationsで

`φ=0, x, x², x²+y, xy`

の五型を尽くす。
rは全H0Ω1(3)源であり、chart上

```
r=P dx+Q dy,
P=A+Bx+Cy+Dxy+Ey²,
Q=F+Gx+Hy−Dx²−Exy.
```

これは任意polynomial degree≤2 one-formの最高degreeのEuler contraction零という
projective global条件そのもので、八coefficientsを落としていない。
flatnessは一つの全polynomial identity

`G−C−3Dx−3Ey+φ_x Q−φ_y P=0`

となる。各φで全coefficientを比較する。

- φ零：D=E零、G=C。r=dψ、ψ全degree≤2 polynomial modulo constant。
  全ψは0,x,x²,x²+y,xyの五型。
- φ=x：D=E=G=H零、F=C。
  r=(A+Bx+Cy)dx+Cdy。
  normal BorelでAを消せ、φを保つy shearでC非零ならBを消せる。
  normal ratio/y scaleで全rは0、x dx、y dx+dyの三型。
- φ=x²：r=(A+Bx)dx。
  normal BorelでBを消し、normal ratioで全rは0又はdxの二型。
- φ=x²+y：r=F dφ。normal Borelで全rを零へ。
- φ=xy：r=C dφ。normal Borelで全rを零へ。

従って全first connectionsは五＋三＋二＋一＋一＝十二fixedへ入る。
全normalizerは前段のpotentialを保ち、negative shear/diagonal inverseで戻せる。
例えばφ=xではy'=sy+tx+uによりC/s、B−Ct/sへ動き、normal Borelで残constantを消す。
φ零ではnormal Borelがrを動かさないため、ψの五型全てを別に残した。

固定planeの安全なsource Bは、Rline二＋πS amplitude一＋constant root一を足した四に、
φ stratum dimensionと元r stratum dimensionを加える。
φ零の五型はB=4,6,7,8,9。
φ=xの三型はB=7,8,9、φ=x²の二型はB=8,9、
φ=x²+yはB=9、φ=xyはB=10。
normal Borelでcanonical r零になっても、その前のF*dφ一parameterをsource側に必ず戻す。
canonical点のstabilizerから全source次元を推測していない。

## 3. R=2L+Mの全remaining first-connection入口

f=a²b、h=a^e b^(2−e)、e=0,1,2。
全pは

`p=(dlogf)/3+αdlog(a/b)+γd(b/a)+δd(c/a)`。

actual c'=sc+tb+uaとtensor scaleで、δ非零はδ1γ0、δ零γ非零はγ1、両零は別に残す。
αは全Cを保つ。従ってe三×exact三の九whole diagonal-connection linesとなる。
q=2p−dlog(f/h)で一意に決まり、rには元八dim H0Ω1(3)全sourceを戻す。
その全linear operator `dr+(q−p)∧r=0` はQQ(α)上有限で、
maximal minorsの全α例外を保持すれば全r-spaceを得られる。
この段階でrの次元又は全normalizer orbitを仮定しない。
normal Borel image(p−q)も、scalar diagonalでそれが零になる全αを別に戻す必要がある。

次はこの全九lineのcovariant-r八次元linear sourceを有限に解き、
actual residual base/normalizerで全rを分類してからglobalK/full30/9へ移る。
rの標本値を選んで全supportを主張しない。

## 4. 射程

R3L十二型及びR2L九line＋全r八dimは全generic triple first connectionsを含む。
raw global quadratic liftのK、P0全Jacobi、higher source及びG非優越は未証明。
特に既知二十三次元componentにrank二plane generic triple-leading実モデルがあるため、
この全branchを局所multiplicityだけで空とは推論しない。
他generic21 degree flags及びm二、reducedRは別の担当sourceである。
