# h5 残実軌道：同一点全周囲の非閉鎖と zero-log 親の閉条件による排除

2026-11-08。PRIMARY `E_multiple_hyperplane`。全 h5 源の既知閉和＋一実軌道への被覆は
`UNIMODULAR_H5_WHOLE_ONE_ORBIT_SEPARATE_REVIEW.md` の独立読解を入力とする。
本稿は残一軌道の親をまだ同定しない。同一点の追加 profile／jet は行わない。

## 1. 同じ forced point の正確な全周囲計算

順序 (xy,xz,xw,yz,yw,zw) で

```
P*=(x,-x²,xz-1,xy-1,-yz+w-y,-xw+xy+y+z).
```

これは旧全 source の n0=n3=a0=0,n1=n2=1 をそのまま戻した点であり、
新しく generic 性を期待して選んだ点ではない。全 n1 n2 unit の正逆 normalizer と
全 translation／最高 GL4／scale の復元は既前稿に固定済み。

`CLASSIFICATION_UNIMODULAR_POLE_SINGLE_ORBIT_ONE_AMBIENT_20261108.json` は元 source から
全六式、全150 homogeneous 係数、strict10／元 projective350、div0、
全五 Pfaffians／全126 quintic 係数を再構成した。実退化因子は h^5。
変更していない既 generalambient engine で全 D350×150、Euler150×25、
Euler＋GL＋scale150×51 の全 entry／完全 kernel／free identity／非零 rank minor を保存し、
元全350行で D M=0 を確認した。

| 行列 | 正確 rank | kernel 次元 |
|---|---:|---:|
| D |101|49|
| Euler |24|1|
| M |47|4|

従って projective tangent 次元は125−101=24、actual orbit 次元は47−25=22。
差二は未解決の法線接方向であり、被約次元、局所枝数、成分性、smoothness の証明ではない。
全 job 1.170600834 秒、最大 RSS 51,948 KiB、180秒／3GiB の契約内で完了した。
**この局所経路は非閉鎖として停止する。別点、追加 minor、higher jet を足さない。**
局所接上界24により、この点を含む全成分の次元は24以下である。
従って証明済みの25以上の成分は親ではない。この推論は新成分性を与えない。

## 2. 固定 h5 chart の affine jet 閉条件

固定 chart の任意の quadratic bivector を P=C+L+Q と書き、
Q の方向微分の linear map を

```
T_Q : V4 -> V4* tensor exterior²V4,  v -> dQ(v)
```

とする。linear part L を一列追加した24×5行列について

```
rank[T_Q | L] <=4                                             (1)
```

は全5角 minor の消滅という閉条件である。

これは affine frame の選び方によらない。linear frame 変更は定義域／値域の可逆変換であり、
tensor scale も非零 scalar にすぎない。translation x -> x+a は
Q を変えず L を L+T_Q(a) へ移すので追加列に既存四列の和を加えるだけである。
特に適切な affine 原点で L=0 となる源の全 affine 飽和とその閉包に (1) が成立する。
ここで「L が T_Q の像にある」自体の閉性は主張せず、弱い determinantal 条件だけを使う。
rank T_Q が落ちる境界も (1) に含まれる。

rank四 open 上、actualD が h^5 の形なら [h] は一意である。
Veronese embedding [h]->[h^5] はその像への同型だから、[h] は actualD から正則に回収できる。
双対射影空間の通常 affine chart で、[h] を固定座標へ送る regular PGL frame を選べる。
その選択の差は h を保つ affine group なので (1) は局所 chart 間で一致し、
h5 rank四 locus 上の相対閉条件を定める。

旧 `legacy/v06/zero_log_component.tex`, Theorem zero_log の定義は

```
C_5H = closure PGL5 * {[h² J+Q(a,b)]}.
```

その全源は actualD=h^5、固定 chart で L=0 である。
従ってその閉包の全 rank四 h5 点にも (1) が成立する。
閉包点近傍で h を正則整列する上記手順があるため、発散する projective frame を
有界 affine translation と仮定する穴はない。

P* の Q と L は

```
Q=(0,-x²,xz,xy,-yz,-xw+xy),
L=(x,0,0,0,w-y,y+z).
```

v=(a,b,c,d) に対する T_Q(v) は

```
(0,-2ax,az+cx,ay+bx,-bz-cy,a(-w+y)+(b-d)x).
```

行 (xz:x),(xw:x),(yz:x),(zw:x) の4×4 minor は −2。
行 (xy:x) を追加列 L とともに加えた5×5 minor も −2 である。
従って rank T_Q=4、rank[T_Q|L]=5。

**P* は既知 R4_ZERO_LOG_5H 成分の全閉包に属さない。**
この結論は genericD≠h5 の他の既知成分を排除しない。
特に beta0 の「有限個の既知親の和」全体の閉包に対する主張へ拡張しない。

## 3. rank-zero 被約台の直接式

Q=X∧Y、X=(x,-y,0,y-w)、Y=(0,0,-x,z) と書く。
h≠0 は Pf1 なので projective rank-zero 点はない。
h0 では actual zero 条件は E4∧X∧Y=0 であり、その四係数は

```
2x²y, -2xyz, x(xy-2xw-z²), y(xy+z²).
```

x非零なら y=0,2xw+z²=0。
x零なら yz²=0 なので y零またはz零。
従って全被約 zero 台は

```
C={y=0,2xw+z²=0}, ell1={x=y=0}, ell2={x=z=0}
```

の和である。全て p=[x:y:z:w]=[0:0:0:1] で交わり、C は ell1 にそこで接する。
従って zero locus は正確に一次元である。
一般元が二次元以上の projective zero locus を持つ任意の既知成分は、
projective zero incidence の properness と fiber 次元上半連続性により全閉包でもその性質を持つので、
そのような成分をこの点の親にすることは出来ない。
本稿は名前付き各成分についてこの前提を新たに証明したわけではない。

## 4. 残る義務

既知閉親への actual polynomial curve／全係数同定が未解決である。
P* が polynomial Darboux coordinates を持つことや全 affine Hamiltonian symmetry は
非線形同型を与える場合があり、PGL 帰属と同一視しない。
全 h5 source の被覆は「既知有限閉和＋この一軌道」までであり、ここを完了扱いしない。
