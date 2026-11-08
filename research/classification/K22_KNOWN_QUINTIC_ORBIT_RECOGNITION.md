# K22 全一軌道の正の認識：既知 nil22-long 既約五次23成分

2026-11-08、PRIMARY G/E。全十括弧を印刷式から手比較した純証明。
CAS、新点、minor、normal jet、profile を用いない。新しい成分は零。
canonical K22 tensor は既知 `R4_NIL22_LONG_IRREDUCIBLE_QUINTIC_23` の稠密23軌道に実 PGL/scale で一致する。
全K22をこの一軌道へ送る作者の whole-normalizer theorem と合わせれば、全K22の既知帰属と像次元23を得る。
whole-normalizer theorem の独立受理はこの十恒等式とは別の入力である。

## 1. 二つの正確な actual tensor

canonical K22 は座標 `(c,a,b,u,v)`、

```
P_K=Q+∂c∧W,
Qab=-av+bu, Qav=Qbu=uv, Qau=Qbv=Quv=0,
Wa=-2a²+4ab, Wb=4ab-2b²,
Wu=-2au+bu-av, Wv=-bu+av-2bv.
```

actual Euler quintic は

```
F*=c*u²*v²+a²*v³+3abuv(u+v)+b²*u³.
```

既知成分の actual dense-orbit tensor P_N は、immutable proof note
`HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_CURVE_RECOGNITION.md` §「第五全源」
に印刷された次の九括弧である。座標 `(H,X,Y,Z,W)`、省略した HX は零。

```
{H,Y}=-6H², {H,Z}=3HX, {H,W}=-3HY+162HZ-45XY,
{X,Y}=9HX, {X,Z}=-2X², {X,W}=-5HZ+2XY-243XZ,
{Y,Z}=-9HZ+3XY, {Y,W}=6HW-3Y²+567YZ,
{Z,W}=2XW+7YZ-243Z².
```

その全ambient D102/G48、smooth unique projective23、および dense orbit23 は既受理 producer/別作者 review に依存する。
ここで新しく ambient 行列を計算しない。

## 2. 非零定数だけを使う actual PGL

target の `(c,a,b,u,v)` から既知座標へ

```
H=27u, X=v, Y=81a, Z=b, W=27c                              (1)
```

と送る。det=27²*81=59049≠0、逆は
`c=W/27,a=Y/81,b=Z,u=H/27,v=X`。
parameter も closed complement も無く、actual linear change である。

E を target Euler、`L=(a+b)∂c-2u∂a-2v∂b` とする。
(1) で P_N を pullback した全十上三角括弧を次表の中央列へ得る。

| bracket | P_N の actual pullback | `(E∧L)` の同じ entry |
|---|---|---|
| ca | `-2cu+9a²-21ab` | `-2cu-a²-ab` |
| cb | `-2cv-21ab+9b²` | `-2cv-ab-b²` |
| cu | `9au-6bu+5av` | `-au-bu` |
| cv | `5bu-6av+9bv` | `-av-bv` |
| ab | `3av-3bu` | `-2av+2bu` |
| au | `2u²` | `2u²` |
| av | `-3uv` | `2uv` |
| bu | `-3uv` | `2uv` |
| bv | `2v²` | `2v²` |
| uv | `0` | `0` |

各行の中央列から右列を引くと、厳密に -5 倍の P_K entry になる。
従って全150係数について

```
(1)^*P_N=-5*P_K+E∧L.                                      (2)
```

これは Euler ambiguity だけの差なので同じ projective Poisson tensor class。
逆 (1) と非零 scale -1/5 で P_K は P_N の同じ actual orbit に入る。
source の quotient-normalization だけから推論せず、全十 entry で actual equality を示した。

## 3. Pfaffian の別の整合性

既知 note の primitive polynomial は

```
F_N=H³Z²+H²WX²+H²XYZ+27HX²YZ+3X³Y²,
actual D(P_N)=75*F_N.
```

(1) を代入すると `F_N=19683*F*`。
actual tensor pullback では determinant を除くから

```
D((1)^*P_N)=75*19683/59049*F*=25*F*.
```

一方 (2) の Euler 項は quintic を変えず、tensor scale -5 の二乗も25。
この actual factor は十括弧の equality と一致する。
primitive witness 方程式を actual Euler Pfaffian の正確な係数と取り違えない。

## 4. Whole source と正確な射程

`PLANE_CONCURRENT_M3_K22_ORE_NONDOMINANCE.md` §§2–7 の whole module theorem は、
全 higher W、A/A'非零の reduced concurrent cubic、全 B/B'/ell を actual inverse shear/scale で除き、
全 K22 をこの一つの canonical P_K orbit へ送る。
その全source theorem が独立受理されると、(2) により全 K22 は既知23成分に入り、
canonical/source像自体は既知 dense orbit23、その閉包は当該23成分である。
MPS40 に属すると呼ばない。これは既知 non-MPS 成分への正の認識である。

E の plane1+1+3残三源や stable二源は、いずれも全source像次元23かつ全 D の本質変数数≤3。
既知 K22 成分の generic F* は本質五変数なので、その閉じた essential≤3 locus に入らない。
従ってそれら23次元の全既約 source 閉包を K22の23成分へ丸ごと含めることはできない。
もし含まれれば同じ次元の既約閉包が全成分に等しくなり、閉じた本質変数条件に反する。
特別な磁気値の低次元 locus が boundary に入る可能性は否定しない。
したがって K22 の線退化は有用な正の全K22帰属を与えたが、残五源の whole parent 探索には24以上の親か別の成分性証明を要する。
