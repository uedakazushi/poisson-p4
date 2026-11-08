# 二つの weak-wall 全源閉包は smooth 二十三次元全成分

2026-11-08。PRIMARY `E_reduced_nonnormal`。producer exact fullambient 算術と純幾何の証明。別作者全行列・全核・minor 再構成は進行中。既126件との相違と台帳加算は統合担当の別義務である。

## 1. Canonical tensor の完全な入力と旧 W0 表示の明確化

`PLANE_M3_TRIPLE_WEAK_WALLS_TWO_ACTUAL_ORBITS.md` は W=W0+全十八の affine solution という意味で W0 を表示した。**その W0 単独は canonical particular ではない**。誤って W0 単独から Ew5 の divisor を再構成してはならない。旧固定 bytes を変更せず、ここで入力を明記する。

座標(a,b,c,z,w)、Π=Qτ+∂c∧W、Qτ は

```
Qab=−(3τ+2)az+w², Qaz=0, Qaw=−3z²,
Qbz=(6τ+1)z², Qbw=9τzw, Qzw=0.
```

実際の完全 particular は次の二つである。

```
τ=4/21: W=(a²+7bw, −2ab, 4az/7+w²/2, aw/7−3bz),
τ=1/12: W=(a²−20bw, −2ab, az/4−w², −aw/2−3bz).
```

十八 module report の `particular_solution`、先の faithful S61 report、今回の ambient report は全てこの完全な入力を使用している。旧 note の D0 table、全 W六 image、実群 transitivity の証明はこの particular を用いており、その数式を撤回するものではない。上記表示の役割を明確にする typed clarification である。

全 W六 source は実 unipotent 六群の一軌道であるという既述の全逆を保つ。今回 integral scale はそれぞれ14、4のみ。新しい任意の source point、slice、specialization は選ばない。

## 2. 全 ambient 計算と成分性

全 quadratic bivector の150係数を、pair順と degree二 exponent の降順で並べる。元 projective Jacobi `E∧[Π,Π]=0` は全五 quad と degree四 monomial の350式を保持する。E を radial field とし、Euler translations E∧V の全25列（rank24）、実 GL5 action の全25列、tensor scale の一列を加えた M51 を用いる。

両 particular に対して、全350×150微分 D、Euler150×25、M150×51を正確に再構成した。結果は両方とも

```
rank_QQ D=102, rank_QQ Euler=24, rank_QQ M=48.
```

全 D*M=0、全四十八次元 D-kernel の元350行への代入零、free-column identity、非零102/24/48 minor を保存した。全行列と全 coefficient/row order は report にあり、selected row だけから上界を推測していない。

raw Poisson scheme の接空間は150−102=48次元である。Euler translation と GL5/scale の真の48次元作用軌道が同じ点を通るため、局所次元は少なくとも48である。従って局所次元と接次元が一致し、その点は scheme として regular、reduced、唯一の全 ambient component 上の smooth point である。Euler quotient は24次元線形空間を全点で戻す product、非零 tensor scale は Gm なので、射影 Poisson scheme の局所次元は48−24−1=23。actual PGL5 orbit の次元も48−24−1=23である。

actual orbit は成分中で同次元の局所閉部分、従って dense open であり、その閉包が二十三次元全 ambient component になる。connected PGL5 は各成分を保つ。全 W六一軌道定理の全逆を戻すと、各 complete weak-wall source 全点がこの成分に入る。normal quotient は零であり、高次障害、jet、C*slice は不要である。

二十三次元という数だけを用いて既成分との相違を推論しない。以前に証明した T39/T40 全閉 graph 除外は引き続き有効であり、他の既成分との比較は別に行う。

## 3. 二成分は互いに異なる：内在的な multiplicity jump

各 actual divisor を

`D=K ab z³+L a w³z+M b w²z²−K c z⁴+E w⁵`

とする。第一は(K,L,M,E)=(45/7,0,15/2,1/2)、第二は(9/2,9/4,−15,−1)。各 D は c-linear primitive、従って C 上絶対既約であり、reduced singular support は唯一の plane S={z=w=0} である。

S の一般点で ordinary multiplicity は三、S 上の multiplicity≥四 locus は正確に La={a=0} と Lb={b=0} の二直線である。これは projective hypersurface から内在的に定まる。一般 La 点では b0≠0 として transverse variables(a,z,w)の四次 tangent cone は

`z²*(K b0*a*z+M b0*w²−K c0*z²)`。

括弧の二次形式は rank三（K,M,b0非零）であり、全 C 上 irreducible、したがって型は二重 plane と irreducible quadric である。

一般 Lb 点では a0≠0 として transverse variables(b,z,w)を用いる。第一 divisor (L=0) の四次 tangent cone は

`z³*(K a0*b−K c0*z)`

であり、三重 plane と別の simple plane。第二 divisor (L≠0) では

`z*(K a0*b*z²+L a0*w³−K c0*z³)`。

括弧の cubic は b-linear primitive（coefficient z² と w³−const*z³ が互いに素）なので絶対既約。simple plane と irreducible cubic の型である。

従って二本の内在的な multiplicity-jump line の一般 tangent-cone factor 型の unordered pair は

```
第一: {2H+irreducible Q, 3H+H},
第二: {2H+irreducible Q, H+irreducible C3}.
```

相異なるため両 divisor は projectively equivalent でない。両 actual tensor は別軌道である。もし二つの全 component が同じなら、その中の二十三次元 actual orbit はどちらも dense open となり、相異なる二つの disjoint dense opens が存在する矛盾。ゆえに両 component は互いに異なる。

## 4. 有限実行と独立確認

事前 gate PASS、資源上限180秒/3GiB。二つの全源を一度だけ計算し3.332秒/RSS149372KiBで完了した。normalized torus diag(s³,s²,1,s⁵,s⁴) と tensor scale s³ は両点を保ち、全150/350 weights と完全核により normal quotient が零であることも同じ行列から確認した。新 rank profile を追加していない。

旧五十一 FunctionDef AST と現40契約 bytes は ambient preappend archive に固定。completed G/E 契約の hash を新 mode へ自動更新していない。raw matrix producer と独立 alternate Fraction reviewer を分離し、別作者の PASS receipt 後にのみ受理済み結果へ統合する。全 prior126との相違・名前・加算は本稿では未完である。
