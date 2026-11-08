# 二つの triple weak wall：全源・全 ambient・相互／prior126 相違の独立 review

2026-11-08、PRIMARY `E_reduced_nonnormal`。別作者の **Fraction 全算術と全幾何の独立 PASS**。
τ=4/21 と1/12の全 weak-wall source は、それぞれ smooth 二十三次元 actual orbit の閉包である
全 Poisson 成分に入る。両成分は相異なり、固定した prior126 の全成分とも相異なる。
新規二成分の数学的採録条件を受理する。catalogue、主稿、Git の変更は統合担当だけが行う。
他の corrected-triple 直接弱上界源や、全 E／全研究の完了を主張しない。

## 1. 実 tensor と typed particular clarification

座標(a,b,c,z,w)、Π=Qτ+∂c∧W。
Qab=−(3τ+2)az+w²、Qaw=−3z²、Qbz=(6τ+1)z²、Qbw=9τzw、Qaz=Qzw=0。
**全 affine module の実 particular** は

```
τ=4/21: W=(a²+7bw, −2ab, 4az/7+w²/2, aw/7−3bz),
τ=1/12: W=(a²−20bw, −2ab, az/4−w², −aw/2−3bz).
```

旧 `TWO_ACTUAL_ORBITS` note の W0 は affine system を解く前の base 部分であり、
それ単独を particular と読んではならない。
この曖昧さを独立 review で指摘し、producer は旧 bytes を保存したまま
`PLANE_TRIPLE_WEAK_WALL_FULL_COMPONENTS.md` §1 に typed clarification を置いた。
全既存 D0 table、W六 image と今回の full ambient は上記の完全な particular を使う。
独立実装も Q と上記四 W entry から元の十括弧を構成し、
integer scale14、4の全公開 actual coefficients と一致させた。

両 canonical quintic は非零 scalar を除き

`D=K ab z³+L a w³z+M b w²z²−K c z⁴+E w⁵`

であり、(K,L,M,E) は第一(45/7,0,15/2,1/2)、第二(9/2,9/4,−15,−1)。
c の coefficient −Kz⁴ と残項の gcd は Ew⁵ のため一。
primitive な c-linear polynomial なので C 上絶対既約、actual tensor の一般 rank は四である。

## 2. 全 W六から一 actual orbit への全逆を再確認

先行の独立 module/factor review は全十八未知 W の affine equation、完全六 kernel、
actualD の全 affine image を元行から確認している。
各 wall の Euler gauge Wz[z²]=Ww[w²]=0 を固定した全 image は

`D0+R(awz³+κw³z²)+A0az⁴+B0bz⁴+C0z⁵+C1wz⁴+C2w²z³`。

κ は7/6と−4/3。
六方向は awz³、az⁴、bz⁴、z⁵、wz⁴、w²z³ の coefficient で独立だから、
同 gauge の全 W六→D map は injective。
chosen W point だけの normal symbol の一致ではない。

η=1−3τ とし、全六 parameter の actual shear

`a'=a+A z, b'=b+C z+ηA w/3, c'=c+ell(a,b,z,w), z'=z,w'=w`

を使う。全 determinant は一で、全 inverse は三角 substitution。
Qab に生じる z² coefficient は ηA−3(ηA/3)=0、他の Q entry も literal に不変。
c-shear は Q core を変えない。
必要なら positive Euler field `(linear(z,w))*∂c` で上の二 coefficient gauge を唯一に戻せる。
これは actual projective action と Euler redundancy のみであり、actual D は変えない。

inverse substitution の D0 に対する awz³ coefficient は−KηA/3、
w³z² は−LA−MηA/3。その比は κ=(L+Mη/3)/(Kη/3)。
w²z³ coefficient は−MC。
両壁で KηM が非零定数なので、任意 R,C2 に対して A,C を全 parameter 上で一意に戻せる。
残る az⁴,bz⁴,z⁵,wz⁴ の四 coefficient は、既知 A,C の多項式を引いた後、
ell の四係数が各−K倍で独立に動かす。
従って全 affine D image は一 orbit、全 inverse は global polynomial である。

群作用は full Jacobi source を保ち、D map は固定 gauge で injective なので、
同 D0 へ戻した tensor は complete particular と Euler を除いて一致する。
全 W六 source は一 actual projective orbit に入る。
前段の全 normalizer inverse を戻しても同じであり、generic dense subset だけの結論ではない。

## 3. 全 ambient D／Euler／実群作用を別実装で再構成

独立 checker `classification_triple_weak_wall_ambient_review.py` は、既済の別作者
`classification_degenerate_quadric_q4_nil31_I_review.py` の standard-library Fraction 疎多項式だけを再使用する。
producer ambient engine、current runner、SymPy、native CAS を import・実行しない。
active finite gate PASS 後、一回の二源 batch で **0.45619秒／RSS21,472KiB**。
7619 shared reviewer と全旧 source/code は不変である。

全十 pair ×十五 degree二 monomial の150座標、全五 quad ×七十 degree四 monomial の350行を一致確認した。
元全350 Jacobi を零と確認後、全150 basis variation から D350×150を作り直した。
Euler translations 全25、GL5 actual action 全25、tensor scale 一列も定義から別に作り、
D・Euler・M51 の全 entry が保存値と完全一致した。

両壁で厳密に

`rank D=102, rank Euler=24, rank M51=48`。

各 matrix の保存 maximal minor は Fraction determinant で非零を再確認。
各 **完全 kernel** は全元行への積零、free-coordinate identity、minor の下界により完全性を証明した。
特に D の全48方向と全350行を残し、D*M51=0も全成分で確認した。
部分 matrix rank 又は少数 coefficient から接空間を推測していない。

全 signed Euler Pfaffian、normalized torus weights(3,2,0,5,4)・tensor character+3、
全150 point coefficients、350 equation weights、D の全 support、全 kernel／group weights も再構成した。
全48 kernel と全48群像の weight multiplicities は一致し、normal quotient は零。
高次障害、normal、jet 又は新しい rank profile は不要である。

## 4. 全成分性と smooth unique 点

raw Jacobi scheme の接次元は150−102=48。
Euler additive action／実 GL5／scalar の真の48次元作用軌道が通るので局所次元は少なくとも48。
局所次元と接次元が一致し、点の local ring は regular、reduced、唯一の full component 上にある。
Euler quotient は全24線形 directions の product、非零 tensor scalar は Gm なので、
projective local dimension と actual PGL5 orbit dimension はともに48−24−1=23。

actual orbit は局所閉かつ成分と同次元、従って dense open。
その閉包が二十三次元全 ambient component である。
§2の全 source inverse により各 weak-wall 全源がこの成分へ入る。
これは必要supersourceの次元上界や同一jet fiberを成分性としたものではない。

## 5. 互いの相違：唯一の singular plane と内在的な二本の線

∂cD=−Kz⁴なので SingD は z=0 に入り、そこで ∂wD=5Ew⁴ から w=0。
逆に z=w=0 では全 partial が零だから、reduced singular support は唯一の plane S。
Sの一般点の ordinary multiplicity は三、leading normal term は K ab z³。
従って S内の multiplicity≥四 locus は正確に二本の直線 La=(a=0)、Lb=(b=0)。
これは hypersurface の内在的 marked locus であり、座標を恣意的に固定した比較ではない。

一般 La 点 b0非零での四次 tangent cone は

`z²*(K b0*a*z+M b0*w²−K c0*z²)`。

括弧は rank三二次形式、従って C 上既約。型は2H＋既約Q。
一般 Lb 点 a0非零で、第一の cone は

`z³*(K a0*b−K c0*z)`、型3H＋H。

第二は

`z*(K a0*b*z²+L a0*w³−K c0*z³)`、型H＋既約C3。

後者は b-linear primitive、coefficient z² と L a0w³−Kc0z³ は coprime なので絶対既約。
両者の二 jump-line 一般 cone 型の unordered pair は相異なり、divisor は PGL同値でない。
両 tensor の23軌道が同じ component なら二つの dense open orbit が交わり同じ orbit になる。
これは矛盾なので二 component は互いに異なる。

## 6. 固定した prior126 の全比較

immutable projection は
`inputs/plane_double_line_m3_triple_two_components_prior126_20261108.json`。
基準 catalogue は commit `9468048f7a9daa051fa57054858aa8f07a1dd61d`、
blob `cd2e6c7bbd66c178777c85209253c5ff4fc10526`。
その全126 IDと全投影fieldを Git blob に照合し、全理由を独立 metadata receipt に保存した。
既成分の成分性・既因子 theorem を理由なく再実行せず、以下の全 inclusion／相違論理を独立確認した。

- rank≤2の六成分：閉 Pfaffian-zero locus にあり、今回の非零 actual quintic を含まない。
- rank四・次元>23の二十一成分：今回の smooth unique projective local dimension23に矛盾する。
- rank四・次元<23の十八成分：PGL-stable closed component が23軌道を含むことは出来ない。
- MPS二十九23成分：既受理 `CLASSIFICATION_MPS23_DIVISOR_FACTORS_20261107.json` の全一般族に proper factor がある。
- nonMPS四十二23成分：既受理一般 divisor は可約又は非被約。
- 既約 essential-three quintic の八23成分：全閉包で essential-variable-count≤3。
- 残る二23成分：Ga9-distinct は二 singular planes、nil22-long は三 singular planes の dense actual23 orbit。

proper-factor locus は正次数の projective multiplication map の proper image の有限和なので閉。
従って一般族での因子 theorem は component 全体に移り、今回の絶対既約 divisor と両立しない。
essential≤3 は全五 partial の coefficient matrix rank≤3という閉 determinantal 条件である。
今回の五 partial の独立性は、順に ab z²、w⁴、b z³、a z³、z⁴ の係数から従い essential五。
八成分も含まない。

最後の二成分の plane 数を閉条件だとは扱わない。
もし今回の23軌道がその既成分に入れば、既成分の既存23 dense orbit と同じ dense open orbitになる必要がある。
PGL同値なら reduced singular support も同値だが、二／三平面と一平面は異なり矛盾する。
既 Ga9 と nil22-long の dense23 orbit・全 singular support はそれぞれの受理済み独立 review を入力とした。
以上は126 IDを全て尽くし、未比較の同次元非平面既約quinticを残していない。

## 7. 正式結論と証明の境界

全元算術、whole W六／全 normalizer inverse、実 orbit lower bound、全接空間、
相互 intrinsic divisor 不変量と prior126 全比較を合わせ、**二つの新しい二十三次元全成分**を独立受理する。
全弱境界源の pointwise E はこれらの成分へ閉じる。
全 Gで次元のみだった他の corrected-triple sourceを、二壁だけで尽くしたとはしない。
新 ID、catalogue 加算、main／global／Git は統合担当がこの frozen proof を採用した後に行う。

算術 phase report、prior126 metadata receipt と本数学証明を、一つの最終 component receipt／manifest に固定した。
旧 false whole-plus coverage、unsolved W0 display、前段の二 weak-wall 存在未主張という当時の限定は
歴史 bytes のまま保ち、今回の成功 proof へ黙って置換していない。
