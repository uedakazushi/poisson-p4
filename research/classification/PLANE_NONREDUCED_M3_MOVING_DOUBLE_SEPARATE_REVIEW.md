# 非被約 R・m三 moving-double 全 degrees(1,0) の別作者検証

2026-11-08。producer と別作者が、全 first connection の正常形・全逆写像を読解し、
元の全係数算術を Fraction／QQ(τ)／有限 Laurent 多項式演算で独立再構成した。
CAS と producer の実行・import は用いない。
本 receipt は **この全 stratum の G 像次元が二十三以下**という命題を受理する。
二つの弱上界源の E 帰属、他 primitive degree、全分類の完了は主張しない。

## 1. 凍結入力と検証範囲

実入力は次の二つの immutable report と五つの数学 note である。

- `CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_20261107.json`
  （SHA256 `56f68485a0da21ee889397ecaeb1f73714e5095650cf9c6c293997600363bd7f`）。
- `CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_20261107.json`
  （`fe69cf2b9c28663427d4aeb402336aa28830677e449960f033041efea572cd0d`）。
- `PLANE_NONREDUCED_M3_MOVING_DOUBLE_FIRST_CONNECTIONS.md`、`PARTIAL_BOUND.md`、
  `ORE_REDUCTION.md`、`ORE_LEADING_ERRATUM.md`、`NONDOMINANCE.md`（同 prefix）。

全 input SHA と独立 checker bytes は
`CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_SEPARATE_REVIEW_20261108.json`
に保存した。独立実行は **11.219496 秒、RSS 53,576 KiB、全 PASS**。
report の producer_bytes は hash 確認だけを行い、実行しなかった。
誤 sign の旧草稿と二つの runtime 停止報告は経緯であり、成功の証明には用いない。
旧 helper df9c は現在の実体と同一とは扱わず、後述 manifest の exact snapshot から復元する。

## 2. 全十一 first connection の被覆

実 D が被約既約 quintic、零 cubic の支持が 2L+M 又は 3L、normal order 三、
leading double row の primitive degree が一、simple row が零、という仮定を保持する。
UFD により simple normal root を z として
`G3=z(Az+Cw)^2`、A,C は独立一次式、C の零線は R の成分、と出来る。
whole normal connection は下三角で、

`q=dlog(f/C²), r=−d(A/C)+(A/C)(q−p)`。

p,q closed が flatness と同値である。p を対数型だけに制限しない。
R=2L+M の全配置は f=a²b、(C,A)=(a,b),(a,c),(b,a),(b,c)、
R=3L は f=a³、(C,A)=(a,b) で尽くされる。
A の C 倍は実 normal Borel で吸収する。

2L+M の全 closed p は
`dlog(f)/3+αdlog(a/b)+γd(b/a)+δd(c/a)`。
閉有理一形式の高位極は exact で、その primitive は一次式/a modulo constants である。
C=a では f r の残る a 極から γ=δ=0、C=b では b 極から α=1/3 が強制される。
3L ではまず `p=dloga+d(P2/a²)`、f r の global 性から a|P2、従って
`p=dloga+γd(b/a)+δd(c/a)` となる。
極係数は primitive A,C の coprimality により零を強制され、chart 内だけの議論ではない。

残る exact 二係数の三場合は、任意の chosen slice への制限ではない。
例えば a'=ta,b'=rb,c'=sc+ub(+va) と、row Az+Cw を保つ実 normal Borel を併用すると
`δ'=tδ/s, γ'=t(γ−δu/s)/r`。
δ≠0、δ=0かつγ≠0、両零の全三場合を、それぞれ (0,1),(1,0),(0,0) へ戻せる。
全 inverse は shear の符号反転と非零 diagonal の逆である。
C=b,A=c は rℓ=sk、C=b,A=a と 3L も対応する二係数の等式だけを課し、
この parameter 作用に必要な自由度は残る。
従って二 whole α line と九 fixed、合計十一の全被覆を受理する。

## 3. 全 global tensor／全30・9／全 parameter 壁

各 β entry から `ωij=3fβij−δij df` を作り、全三成分の polynomial 性と
Euler horizontality を再検証した。`Vij=curl(ωij)/3` と逆式 `ωij=−x×Vij`、
trace-free 性も元係数で確認した。raw mixed bracket は
`{xk,n_i}=−Σj Vij,k n_j`、base block は gradient f である。
全 K=k0 z²+k1 zw+k2 w² の D2 零条件三十係数を再構成し、
rank／augmented rank と原行積により一意 global K を確認した。
全十一 P0 は actual D=0、全三百五十 projective Jacobi=0。

全 P1 三十方向、P2 九方向を初めから戻し、全350×30／350×9 operator、
actual D3 の全24×30 map、`D3∈C*z(Az+Cw)^2` の全二十三 annihilator 行を再構成した。
各 complete kernel は、元の全行との積零、free coordinate の単位行列、
指定された正方 minor の正確な非零 determinant で証明した。
actual six shear も元 tensor から再構成した。
二 whole α line の全 generic minor と kernel 分母の全根を独立に線形因子除算し、
両 line 合計十四の元 fullsource 壁を全て再計算した。
今回全根は有理であり、非有理複素根を無視したものではない。

七 whole 源は全 kernel の D3 像が零で m三対象が空。
R2L,C=a,A=c の whole α は全壁込み像二十三以下、
R3L の exact-C 源は二十二以下である。
上界は whole first-symbol source B と k1,k2 から `B+k1+k2+5`。
固定 α 壁は B を一減し、P2 の非斉次条件は P1 ごとに空又は kerJ2 の affine 部分集合なので、
二次 compatibility を捨てた必要 source で十分な上界となる。
五は plane の Gr(3,5) 六を戻し projective scalar 一を引いた値である。
B の whole parameter count は producer §6 の root data／exact 二係数を全て保持しており、
一般 orbit 次元を特殊点へ外挿していない。

残る二 source の粗上界三十三／三十二を、そのまま完了とはしない。

## 4. 残二 fullsource と三 Ore module の全逆写像

残二全 source の P1 kernel 十三と **P2 全九**を戻す。
最初の Euler gauge は `L=(-∂c P_ac,−∂c P_bc,0,0,0)`。
全 c 依存が消え、次の c 方向 Euler gauge で Wz の z²、Ww の w² を除ける。
全 tensor は projective に `Q(a,b,z,w)+∂c∧W(a,b,z,w)` となる。
Q は u と全 central quadratic F の四係数、W は十八係数である。
元全二十二係数からこの二十二係数への joint linear map は、両源とも
**determinant −1**。全 source の literal Euler gauge、c 非依存、全 Q/W 係数一致を検証した。
従って二十二の自由な係数空間を chosen Q slice から推論していない。
旧 operator 壁で空でないものも kernel dimension 十三を復活している。
generic kernel は polynomial、free coordinate が常に単位行列なので各壁でも独立十三本、
同じ完全 kernel と一致する。この議論で joint map の全 t 被覆を確保する。

R2L では A=2/3−2t、B=A+1 とし、全 Q は producer の五括弧の式となる。
全 strict Q Jacobi は恒等零。**訂正済み** leading は

`D3=(3/2)u A(A−1) z(bz+aw)^2`。

R3L は `D3=6u z(bz+aw)^2`。
従って u 非零、R2L は A(A−1) 非零が必要である。
除いた t=−1/6,1/3 は元 full D3 map でも全零である。
実 diagonal の tensor character sa/sc、sb=sa sw/sz により
`u'=sa²u/(sc sz)`、R3L では sa=sb,sw=sz で同じ正規化が出来る。
全 Q と W0 の pushforward 恒等式を別作者 Laurent 算術で確認した。

実 positive-weight shear a'=a+αz+βw,b'=b+γz+δw の full Q 恒等式から、
R2L central change は `u[(A−1)γz²+A(α+δ)zw+Bβw²]`。
A(A−1)B 非零なら全 F を消せる。
B=0、t=5/6 では F=ηw² を保持し、η零／非零を分ける。
u=1 を保つ残 diagonal は η'=sa η/sw なので非零 η を一へ戻せる。
従って whole t,F=0 と一 extra t=5/6,F=w² で尽くされる。
R3L の central change は
`u[(2γ−δ−α)z²+(2δ+2α−2β)zw+2βw²]`。
α,β,γ の三列 determinant は八、全 F を消せる。
全逆 shear は符号反転で、W は十八全係数へ regauge され、c を移していない。
これで三 Ore module への **全 source** の actual 正規化を受理する。

## 5. 全 affine module と全例外の空性

三 module の full350 affine W 十八 operator を定義から復元した。
W 同士の Jacobi 項は c 非依存な共通 ∂c exterior factor のため恒等零であり、
この affine operator が元全 projective 条件である。
全 matrix と augmented matrix の rank／全 complete kernel／particular solution を元行で再検証し、
actual D の全 affine parameter 係数を再構成した。
proper factor は全 particular と全 kernel 方向に対する正確な多項式除算で確認した。

| module／t | full W 核 | 全 actual D の proper factor 又は矛盾 |
|---|---:|---|
| R2L whole generic | 4 | z |
| t=4/3 | 5 | aw+bz |
| t=5/6,F=0 | 6 | z |
| t=7/12 | 5 | z |
| t=1/12 | 5 | z |
| t=−1/6,1/3 | 元 source | leading 全零 |
| t=5/6,F=w² | 5 | rank13、augmented rank14 で不整合 |
| R3L fixed | 4 | z |

特に t=4/3 の核五を generic 核四や z 因子に置き換えていない。
全 generic minor の零点、kernel／particular／factor の係数分母の根は表の全六 t に支持され、
各根を元 matrices へ戻した。線形因子の繰り返し除算が定数まで落ちることにより、
隠れた複素壁が無いことも独立に検証した。
全て不可約 D の必要条件に反するので、残二 Ore whole source は空である。

従って十一 first connection は空又は全 projective 像二十三以下。
この G 結論と、二弱上界源の未解決 E 帰属を区別する。
本検証は新成分、他 primitive flag、m二／m四、reduced R の結果を追加しない。
