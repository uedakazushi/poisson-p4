# twisted cubic cone・全base軌道とvertex fiberの別作者確認

2026-11-07。producer は generic_exhaustion、reviewer は plane_collision_review。
対象は `CUBIC_CONE_ABELIAN_CONORMAL_NONDOMINANCE.md` の最終 V2 証明。
**判定：全 rank零・既約被約 quintic 枝の PGL 飽和像 dimension≤24 として PASS。**
初回 ORBITS 報告の誤った semisimple ラベル又は universal boolean は使用しない。
cone の全十八変数 Poisson 源そのものの dimension≤9 も主張しない。
base=0 に実在する十二次元 rank二源は、そのまま保持する。

## 1. 全三十源からの必要十八源

twisted cubic cone の三二次 ideal は

`q=(x0x2−x1²,x0x3−x1x2,x1x3−x2²)`。

Hilbert–Burch 列と Euler 全源補題から、rank零の全 bivector は `Σqg Ag` の三十係数である。
scroll の別レビューと同じ nonplane multiplicity二・Frobenius first-jet 議論を再確認した。
必要なのは smooth 一般点での abelian conormal と、それから全 first partial の
polynomial 消滅を得る方向だけである。cone vertex での smoothness を仮定しない。

独立 checker は `(s³,s²t,st²,t³,c)` の dominant parametrization へ引き戻し、
全78×30行と labels、rank12、非零minor、完全30×18kernelを Fraction で再構成した。
元の三十 basis と新十八 basis の双方から全350 Jacobi・全126実五次係数を作り直し、
保存された全係数と一致した。

## 2. 六base・十二vertexと、九式の全支持

全kernelを実 tensor の support から検査すると、一意に

`Π=P(x0,...,x3)+∂c∧V(x0,...,x3)`

と分かれる。P は六次元、V は全 `Σqi Vi` の十二次元である。
全350行は pure base 二次又は base×vertex に分かれ、vertex² の式は無い。
全六括弧表示、350個の base 行、350×12の全 symbolic fiber matrix を再構成した。

九本の base 方程式について、元350行との exact rowspace の一致と、
保存した元行 functional を全て確認した。
その被約支持が次の和で尽くされる手証明も確認した。

```
L3: B=0, F=2A, D=−C,
Q2: C=0, E=6B, F=3A, AD=4B².
```

A≠0 では最初の積から F=2A 又は3A。
前者は B=0,D=−C、後者は E=6B,C=0,AD=4B² を強制する。
A=0 では F²=0、BE=0、12B²−4BE=0 から、被約支持で F=B=0。
残りは `C(C+D)=0,(C+D)E=0`。
C+D=0 なら L3、そうでなければ C=E=0 で Q2 に入る。
逆向きの全九式への代入も exact に確認した。
従って A=0又は nilpotent 境界を、一般開集合の外に捨てていない。

## 3. 実SL2作用と、正しい三つの全軌道stratum

twisted cubic parametrization の SL2 が与える h,e,f の ambient 線形場を用い、
元の六 base basis すべての Lie derivative を独立に微分した。
得られた full6×6行列、L3とspan(Q2)の相補基底、全block行列は V2 と一致した。
両blockは h の固有値2,0,−2と非零の隣接 e,f arrows を持つ。
それぞれ実際の Sym² 表現であり、単に保存された invariant を信頼したものではない。

L3 の invariant は `E²+4AC`、Q側は `AD−4B²`。
各3次元blockで h,e,f が独立な三次元 algebra を張り、非退化二次形式を保つことも確認した。
従って SL2 はこの形式の connected special orthogonal group を与え、source scalar と合わせて
非零 discriminant の開集合と非零 nilcone が、それぞれ一軌道となる。

全非零 base は厳密に次の三 stratum である。

| stratum | base dimension | 正しい代表(A,B,C,D,E,F) |
| --- | ---: | --- |
| L3、invariant非零 | 3 | `(0,0,0,0,1,0)` |
| L3、非零nilcone | 2 | `(0,0,1,−1,0,0)` |
| Q2、非零nilcone | 2 | `(1,0,0,0,0,3)` |

E=1 の invariant は1である。初回報告の A=1,C=E=0 は invariant零であり、
semisimple 代表ではない。
checker はその初回報告の不変SHAと ERRATUM の corrected-report 指定を検証し、
初回の全stratum結論を明示的に無効として扱う。

## 4. 全vertex fiberと次元上界

三代表における全350×12 matrix を元の symbolic matrix から再評価し、
全kernelの代入零、kernel独立性、全matrixのrank、非零rank minorを確認した。

| stratum | fiber rank | 完全kernel | fixed source上界 | PGL像上界 |
| --- | ---: | ---: | ---: | ---: |
| L3 invariant非零 | 6 | 6 | 9 | 24 |
| L3 nilcone | 6 | 6 | 8 | 23 |
| Q2 非零 | 9 | 3 | 5 | 20 |

SL2 と scalar による変換は全fiber線形系の同型を与えるので、rank は各全軌道上で一定。
したがって三つの任意標本からの一般性推論ではなく、全stratumの dimension 上界である。
base=0 は `∂c∧V` と decomposable なので rank≤2。
今回の rank四枝には入らず、別の rank二成分分類を仮定しない。

cone preserving group は binary-cubic PGL2、vertex相対scaling、四つの vertex translation により
忠実な八次元群を含む。全必要 source と上のstrataを保つ。
`PGL5×X→P H0(Λ²T)` の fiber に、この八次元群と scalar の九次元が自由に入る。
個々の tensor の stabilizer が零であることは要求しない。
よって projective image は `24+d−8−1` 以下となり、表の上界が従う。

## 5. 不変保存物と結論の限界

- `code/classification_cubic_cone_abelian_review.py`
- `research/classification/jobs/cubic_cone_abelian_review.json`
- `reports/CLASSIFICATION_CUBIC_CONE_ABELIAN_REVIEW_20261107.json`

算術は既受理の scroll reviewer の Fraction/疎多項式 primitives のみを再利用し、
producer、producer helper、SymPy、CAS は import していない。GB 再実行も無い。
事前 gate 後0.367秒で終了した。producer/source/V1/V2/ERRATUM の bytes は変更していない。
最初の preflight は JSON の文字列 matrix と Fraction を正規化せず比較したため停止し、
その実装上の比較を修正した後に上記一つの最終報告を保存した。

smooth scroll の独立上界23と合わせ、degree三の全非退化 rank零曲面枝は G 一般元から排除される。
これは全 Poisson source の E 帰属、rank四 source の不存在、新成分数、他degree又は
G 全体の完了を意味しない。
