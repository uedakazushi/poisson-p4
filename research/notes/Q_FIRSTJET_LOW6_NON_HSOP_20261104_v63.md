# 低6不変式はHSOPでない：Hilbert係数−3と正次元zero fiber

2026-11-04。標準的なHochster–Roberts定理と、凍結済みの**完全な**character次元を
組み合わせる小さい構造的obstructionである。新しいbeta/source点、source26展開、
不変式環の再構築、GB、定義式のdegree chaseは行わない。

## 1. Ringと六つの実際の不変式

```
G=Sp4,C, W=W16=V_(omega1+omega2),
S=C[W]=Sym(W*), R=S^G.
```

Rはconnected・positively gradedな有限生成domain、dim R=6である。
Sp4は標数0でlinearly reductive、Sはregular polynomial ringなので、
**Hochster–Roberts定理によりRはCohen–Macaulay**。

既存の全係数恒等式とexact六次minorで代数的独立性が証明されている六つを

```
f=(F4,p1,p2,p3,p4,N12),    degrees=(4,8,8,8,8,12),
A=C[f] subset R
```

とする。Aは六変数のweighted polynomial ringである。
五つのpencil係数には `F4²=400p0+1600p4` があり、四つのp8を選ぶ理由はこれである。
Rの次元6やInv8の線形独立性だけから、上の特定の六つの独立性を推定してはいない。

## 2. H26=45は候補subsetの大きさではない

`Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW` は、元の16 weightsの
Sym冪とC2のWeyl denominator全8項を26次まで独立再構成している。
全57候補のrankが45という下界だけでなく、characterの**全空間上界45**と一致し、
Inv26全体のQQ基底が45本であることを受理している。
従ってC上でも同じ完全次元であり、候補不足をHilbert次元と誤認していない。

保存した全characterから必要な値は

|d|2|6|10|14|18|22|26|
|---|---|---|---|---|---|---|---|
|dim R_d|0|0|0|3|9|21|45|

今回のcheckerはSym冪や不変式を再計算せず、保存した全8項の符号付き整数和を
再確認し、独立レビューとその全character file hashを結び付ける。

## 3. Hilbert numeratorに必然的な負係数

`H_R(t)=sum dim(R_d)t^d` とし、候補parameter degreeに対応する多項式を

```
D(t)=(1−t4)*(1−t8)^4*(1−t12)
```

と置く。26次までに関係する係数は

```
D(t)=1−t4−4t8+3t12+7t16−2t20−8t24+O(t28).
```

従って

```
[t26] H_R(t)D(t)
 =H26−H22−4H18+3H14+7H10−2H6−8H2
 =45−21−4*9+3*3
 =−3.                                         (1)
```

これは有限の整数convolutionそのものであり、Molien級数全体の閉形式や未計算の
高次数係数を仮定しない。

もしfがHSOPなら、Cohen–MacaulayなRではfはregular sequenceである。
各regular homogeneous elementの短完全列を順次使えば

```
H_(R/(f))(t)=H_R(t)D(t)
```

となる。左辺はArtinian graded C-algebraのHilbert級数で全係数が非負。
(1)は矛盾。従ってfはHSOPでない。さらに計算はfの個別係数を使わないので、
**Rにはdegree type (4,8,8,8,8,12) のHSOPそのものが存在しない**。

## 4. Generically finiteであるがfiniteでなく、zero fiberは正次元

dim A=dim R=6なので、Frac R/Frac Aは有限拡大であり、商写像

```
pi: Spec R -> Spec A=A6
```

はgenerically finite dominantである。しかしRがA上module-finiteなら、fはHSOPに
なってしまうため不可能。従って **このgenerically finite写像はfiniteではない**。

さらに、positive gradingを使うと

```
RがA上finite iff R/(f)がC上finite-dimensional.
```

右から左は、quotientの有限homogeneous基底を持ち上げ、任意の高次数元を
sum f_i*s_iへ落とし、そのs_iの次数が真に小さいことを使って次数帰納すればよい。
A上の有限性を事前に仮定したNakayamaではない。

従ってR/(f)は無限次元。有限生成graded C-algebraでdimension0ならArtinianかつ有限次元
なので、**pi^-1(0)=Spec R/(f) は正次元**である。

これは商空間Spec Rでのlow6 zero fiberであり、「全ての正次数不変式が零」という
通常のGIT nullconeと同じではない。正次元zero fiberはvertex以外のclosed pointを持つ。
reductive affine quotientの全射性により、low6は零だが別の高不変式は非零となる
Wの点が存在することも従う。ただし今回はその点を一つも構成・評価していない。
zero fiberの正確な次元、成分数、最小の追加invariantは未決である。

## 5. 残余曲線の境界とNEW定義式への含意

low6=0なら規格化恒等式からp0も零なので、p0,...,p4全て零である。
既存のordered regular-semisimple chartでは

```
F4=20(x²+y²), p4=x²y²,
```

だからF4=p4=0はx=y=0を強制し、最初からxyDelta!=0というchartの外にある。
よって今回の正次元fiberは、generic residual curveの種数3やm-degree8の計算で
除かれた**本当のspectral boundary**に存在する。N12=2p0*mのm自体もp0=0では
有理表示の分母を失う。generic coverをそのままoriginへ延長したと扱わない。

特に「generic ambient field degree」と「NEW first-jet hypersurfaceのgeneric spectral degree」
だけをHilbert multiplicityに代入して、NEWの定義式次数を30等へ決定する方法には、
未証明どころか今回否定されたfinite-HSOP仮定が含まれる。
これはdegree30の方程式が存在しないという否定ではない。その推論が成立しないという意味である。
定義式の構成・正確なdegreeは引き続きopen。

次の構造的な対象は、low6 zero fiberに残る高不変式、あるいはこのspectral boundaryを含む
適切なfinite parameter systemとpole/divisor dataである。今回、新しい点・追加HSOP探索・
source26展開を開始したわけではない。

## 6. 標準定理の引用と保存

M. Hochster and J. L. Roberts,
*Rings of invariants of reductive groups acting on regular rings are Cohen-Macaulay*,
Advances in Mathematics **13** (1974), no.2, 115–175.
DOI: <https://doi.org/10.1016/0001-8708(74)90067-X>。
書誌は出版社metadataと[著者Robertsのpublication list](https://www-users.cse.umn.edu/~rober005/Publications.html)
で確認した。ここではchar0のlinearly reductive群・regular polynomial ringへの標準的適用を
用いるもので、原論文の証明を今回再現したという主張や新規性の主張ではない。

専用source `code/q_firstjet_low6_non_hsop_20261104_v63.py`。
同名reportsに全27個の完全Hilbert係数、分母多項式の全係数、convolution、入力hashを保存。
計算は小さい整数和・積だけであり、R rebuild、CAS、点、不変式source評価は行わない。
