# constant TRIPLE：被約cubic上の全affine葉は次数六

2026-11-08、PRIMARY E_multiple_hyperplane、担当 `/root/exceptional_coverage`。
**純手展開。新point、operator、rank、変形jet又はCASの実行無し。**
T37側の全source／proper invariant-leaf条件は別担当の証明であり、本稿はtarget側だけを証明する。

## 1. 全canonicalのrank二restriction

全四parameterを実affine/PGLで一actual orbitへ送る受理済み正常化を使う。
座標(h,t,x,y,r)、h-row零で、canonicalの六括弧は

```
Ptx=xh/3+yh,
Pty=−xy+yh/6−rh/3,
Ptr=−3y²+xr−5yh/4+rh/2,
Pxy=h², Pxr=−3xh+3h²/2, Pyr=3x²+3yh.
```

実D=h² C、C=x³+3hxy−3h²y/2+h²r。Cは絶対既約である。
h=1のC上では

```
r=−x³−3xy+3y/2.
```

これはA³_(t,x,y)との同型を与え、特別点又はparameter unitを除かない。
CはPoisson divisorなので、restrictionの全三括弧は

```
{t,x}=x/3+y=:A,
{t,y}=(x³−y)/3=:B,
{x,y}=1.
```

第二式では−xy+(xy)が相殺し、y/6−y/2=−y/3となる。
{x,y}=1なので、この全affine chartでrankは厳密に二である。

## 2. 核一形式は全域でpolynomial exact

順序(t,x,y)のskew matrixの核は(1,−B,A)。従って

```
alpha=dt+(y−x³)dx/3+(x/3+y)dy=dK,
K=t+xy/3−x⁴/12+y²/2.
```

閉性はdx dyの二項が相殺することで直接確認出来る。
さらにdK/dt=1。全k∈CについてK=kのfiberは

```
t=k−xy/3+x⁴/12−y²/2,
r=−x³−3xy+3y/2
```

というA²_(x,y)のgraphである。従って全finite fiberは滑らか・既約であり、
rank二分布はfiberの接束に等しい。各fiberは一つの連結な複素symplectic leafである。
k=0を含む全finite値を扱い、generic kだけの計算ではない。

Cに含まれ、このchartの一点を通る任意のirreducible invariant surface Wを取る。
その滑らかな一般点では接平面がrank二分布と一致し、d(K|W)=0。
標数零でWは既約なのでK|Wは定数kとなる。
同次元であるからWのこのchart内の閉包は上のfiber全体。
従ってprojective invariant surfaceもこのfiberのprojective closureに一致する。

## 3. 全finite葉のprojective degreeは六

葉のaffine graphのprojective closureをS_k⊂P4とする。これは既約surfaceである。
一般の二hyperplaneを、t,rの係数の2×2 minorが非零の開で解くと、交わるP²は

```
t=A0*h+A1*x+A2*y,
r=B0*h+B1*x+B2*y
```

と書ける。これだけがplane選択の開条件であり、source点やkを変えていない。
P²上でCとK=kのhomogenizationはそれぞれ

```
C3=x³+3hxy+(B2−3/2)h²y+B1h²x+B0h³,
G4=−x⁴+4xyh²+6y²h²+12h³(A0h+A1x+A2y)−12k h⁴.
```

一般planeではproper intersectionである。
h=0ではC3=x³、G4=−x⁴なので、唯一の共通無限遠点はp=[h:x:y]=[0:0:1]。
y=1としbeta=B2−3/2と書く。C3の最低項はh(3x+beta h)で、
二つの接線はbeta=0の場合も異なる。従ってpにおけるC3は二つの滑らかな形式的枝を持つ。

* h=0を接線とする枝では、C3の次数三からh=−x²/3+O(x³)。
  この枝でG4=−x⁴+6x⁴/9+O(x⁵)=−x⁴/3+O(x⁵)となり、交点重複度は四。
* もう一つの枝はhをuniformizerに持ち、x=−beta h/3+O(h²)。
  この枝でG4=6h²+O(h³)なので、交点重複度は二。

よって局所交点数は4+2=6。
二つの非零leading coefficient −1/3と6はA_i,B_i,kに依存せず、
beta=0又はk=0を捨てない。Bézoutからaffine交点長は3*4−6=6となる。

S_kの無限遠境界は高々一次元なので、一般P²はその境界を避ける。
したがって上のaffine交点長はそのままdeg S_kである。
C3/G4の無限遠の共通点をleafの余分な点と数えたり、
complete intersectionの次数十二をleafの次数と取り違えたりしていない。

```
全k∈Cについて deg(S_k)=6.
```

## 4. bounded invariant cycleに対する最弱の帰結

実D=h²CのC∩{h≠0}の一般点pを取る。
Dに含まれるeffective二cycleがpを通り、各irreducible surface componentがPoisson invariantなら、
pを通るcomponentはCに含まれる。
§2–3からそのcomponentの次数は六であり、cycleの全次数も少なくとも六である。
非被約multiplicity又は別componentを加えても次数は減らない。

ゆえに、このcanonicalには
「Dの一般点をdegree≤5のPoisson invariant surface cycleで覆う」
という性質が無い。全四parameter源は一actual orbitなので同じ結論を持つ。
PGLはsurfaceの次数とPoisson不変性を保つ。

別担当がT37全generic sourceにこのdegree≤5被覆を構成し、
proper Hilbert／cycle incidenceによって全rank四境界へ延長するなら、
このtarget計算は全T37非所属を与える。
そのsource完全性・閉性の証明を本稿で済ませたとはしない。
特に正規化や葉の平坦族を仮定するだけで境界へ外挿しない。

T39非所属は別の全graph陰性証明書による。
本稿単独では新成分性、全prior比較又は全E完了を主張しない。
元source、共有code/job、既report、global台帳及びGitは変更していない。
