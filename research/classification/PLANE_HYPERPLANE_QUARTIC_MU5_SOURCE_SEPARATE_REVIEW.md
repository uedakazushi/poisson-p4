# μ5 marked quartic：全配置・全六磁気源・known103境界の別作者検証

2026-11-07、独立 reviewer `/root/quartic_review`。**PASS**。
対象は type6 f0 cotangent、D=HF 被約、F 絶対既約四次、μ(F)=5 の marked branch。
入口 τ(D)=12 より H∩F は一点、全 D の特異点は quasihomogeneous である。
producer の全五 profile の marked list、全五有限源と A3+A2 k=−1 control を査読した。
全周囲の四 componenthood は root の別作者検証
`CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_AMBIENT_REVIEW_20261107.json`
に分離する。**A5 smooth-hyperflex の known-parent membership/成分性は未知**。

## 有限算術の独立再現

checker `classification_marked_quartic_mu5_review.py` は frozen 別作者 MU6 checker の
Fraction・疎多項式演算を利用し、producer/CAS は import しない。
gate PASS 後、全六 case を 0.156 秒、17,108 KiB で再構成した。
実行 bytes は `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_SOURCE_REVIEW_20261107.json` 内にも保存する。

全六で actual H→h GL3及び逆、全28×18 gradient coefficient map、rank16 小行列式、
完全な二次元 kernel、h-Casimir 二 anchor、Saito determinant を照合した。
固定 basis だけでなく**全 GL2 frame 四係数及び磁気 R 全六係数**を独立不定元として、
全十 strict Jacobi（従って全350 projective Jacobi）と実 Pfaffian を全係数確認した。
発散 N は全源で base を殺す fiber 成分のみ、係数は base 一次式。
よって N²=0 は構造的恒等式、canonical rank2/h双対像外も全六で確認した。

全 T は rank5。保存の非零5minor と左 cokernel v を使い、
全磁気を `R=(v m/v_j) mon_j` へ送る**正則な6×6shear inverse**を構成し、
`T B=N_projection−Id` を全六係数で再生した。
未知磁気係数の反転は不要。I≠0 は fiber homothety と tensor scale で canonical へ移し、
I=0 は実多項式族 R=t mon_j の端点で戻す。GL2 frame の変更も shear 定義域を可逆に変えるだけ。
従って各全六磁気 source の PGL 飽和閉包は canonical 一軌道閉包に等しい。

base marked-pair automorphism は有限。
A4+A1、2A2+A1、A3+A2 は normalization 上の三点以上の cusp/node 逆像を保つため有限。
A5 は二逆像固定の r-scale が y section の s³t/st³ 比を保つ条件で有限。
D5 は唯一特異点と異なる二接線を保つ行列を F の h4,y4,h²yw と比較すると base scalar だけ。
intrinsic な fiber kernel と唯一 H を保つ tensor stabilizer は base 恒等に還元でき、
anchor equality が A=σId、I≠0 が σ²=σ、従って σ=1 を強制する。
残りは kerT の一次元。各 canonical orbit は正確に23次元である。
I=0 の全補集合を含む source closure をこの dense orbit closure と同定する。

## 個々の実曲線の presence と型認識

各 degree4 parametrization の全 degree3 implicit coefficient map は rank10、共通 basepoint は無い。
F の pullback は恒等零。よって image は絶対既約四次、map は birational である。
これは必要な全リスト分類を先に仮定しない、個々の actual presence の検証である。

* A5：`(s²t²,st(s²−t²),s4+t4−4st(s²+t²))`。
  二 smooth branch の局所 h(y) を逆級数で解き、差の最初の非零項は三次、従って A5。
  marked H=w+6h の pullback は (s−t)4、像は smooth。
* A4+A1 node-inflexion：`(s4,s²t²,2(s³t−st³)+t4)`。
  cusp の first characteristic odd5、node [1:1:1] は全 gradient零・Hessian determinant非零。
  H=h−w の pullback は (s−t)³(s+t)、node 上の接触3+1。
* A4+A1 smooth-hyperflex：`(s4,s²t²,t4+4st(s²+t²))`。
  cusp odd5、node [1:−1:1] の Hessian非零、H=h+6y+w は nodeを避け、pullback (s+t)4。
* D5 smooth-hyperflex：`(s³t,s²t²,−s4+4st³−t4)`。
  一枝は A2、他枝は smooth で cusp と transverse、従って D5。
  H=4h−6y+w の pullback は −(s−t)4、像は smooth。
* 2A2+A1 smooth-hyperflex：`(s²t²,s³(s+t),t³(s+t/16))`。
  両端 cusp は odd3、node [1:−1/16:−1] は gradient零・Hessian非零。
  H=6h+16y+w は nodeを避け、pullback 16(s+t/4)4、像は smooth。
* A3+A2 k=−1：`(s²t²,st(s−t)²,−s4+4s³t+t4)`。
  tacnode二 smooth branch の局所差は二次、他の cusp は odd3。
  H=h の pullback は s²t²、tacnode common tangent。

各局所 δ の和は3なので、他の特異点や複数 image branch は無い。
µ5四候補の異なる profile、同 A4+A1 の node/H incidence はこの直接認識で証明される。

## 全 marked-family 必要性の読解

有限代表だけから全 source を推測していない。
全 μ5 profile A5、D5、A4+A1、A3+A2、2A2+A1 について、
最初に全 global parameter を保持し、それから binary pullback の一 root 条件を解く。

A5 では二逆像を0,∞に置き、共通接線 h=s²t²、transverse section y=st(s²+t²)。
第三 section の両端係数の一致が branch 接触≥3を表し、target shear 後
`w=s4+t4+k(s³t−st³)`、k≠0 が全 family。k=0 は二対一 mapで integral quartic normalization でない。
H が A5点を通る場合は共通接線だけだが D 三枝の pairwise contact(3,2,2) は weighted homogeneous でない。
H が点を避ける場合は全一 root比較が u4=1、k=2(u³−u) を与える。
u=±1 は k0補集合、u=±i は同一 marked orbit。唯一 A5 smooth-hyperflex 候補を得る。

A4+A1 は cusp逆像0、node二逆像±1を固定した全
`(s4,s²t²,k(s³t−st³)+t4)`、k≠0。
A4 tangent H は intersection4だが、A4 weights(2,5) と同時 homogeneous な smooth branch の intersection は2又は5、
従って τ12条件から除外される。
nodeの唯一交点条件は roots(3,1)、全係数比較で k=±2、符号は実 r→−r で交換。
node/cuspを避ける一 root条件は u4=1,u²=−1,k=∓4i、同じ rational smooth-hyperflex pairへ移る。
これで二 marked orbit の他に kを残さない。

D5 の全 marked family は `F=h²yw+h4+dhy³+y4`。
h4,y4 の元係数は零なら F が可約なので、両方の unit を正規化する。
HがD5を通ると cusp tangent h 又は smooth-branch tangent y の二つだけ。
h枝は全 τ12正規化 d0によって既 A0へ被覆され、MU6別査読で受理した。
y枝は A2 cuspに intersection2を持つ別 smooth branchが二つ必要になり、
A2 weights(2,3) の同時 homogeneous axis条件に反して除外される。
D5点回避の smooth hyperflex は全 binary式
`−Cr4+Ar3+Br²−Cd r−C=−C(r−u)4` から
u4=1,d=−4u³,A=4Cu,B=−6Cu²。四乗根の実作用で一 orbitである。

A3+A2 全 family は `k≠0,1` の
`(s²t²,st(s−t)²,k s4+2(1−k)s³t+t4)`。
k0は basepoint、k1 は tacnode曲率差を失い、対象外として別に保持する。
smooth hyperflex の全係数比較は ku4=1 と
`(u−1)³(u+1)=0` を与え、u1は cusp、u−1は k1なので不可能。
A2 tangent contact3は不足し、唯一交点4となる H は tacnode common tangentだけ。
全 k は下節の known103認識に戻る。

2A2+A1 全 family は
`(s²t²,s³(s+t),t³(s+λt))`、λ≠0,1/4。
λ0は basepoint、λ1/4は node二逆像が衝突して3A2となり、既 MU6枝へ戻る。
nodeの二逆像は u+v=−1,uv=λ。
node唯一交点の(3,1)係数比較はu0又はv0を強制しλ≠0と矛盾。
cusp tangent contact3も不足。一つの smooth root を持つ条件は
u=−1/4,C=B/16,λ=1/16,A=3B/8となり、唯一 hyperflex pairに尽くされる。
以上で全五有限 source、既 A0、全 A3+A2→known103 以外の μ5 marked case を落としていない。

## A3+A2 全 k、特に k=−1 の actual source 閉性

polynomial family を

```
V=hw−y²−4hy+(k−3)h²,
Fk=V²−(k−1)y(y+2h)V+(k−1)²h²y²
```

と置く。全 normalization への代入零と k−1 specialization を独立展開した。
k≠−1 で b=(1−k)/(2(1+k))、a=b² とすると、既知 F3(a) への actual GL3 は

```
H=16b4 h,
Y=4b² y,
W=−6(2b+1)w+24(b+1)y+(12+48b)h.
```

det=−384b6(2b+1)≠0。これは normalization の実 Möbius
`(s0,s1)=(b(s+t),s−t)` から得られる。
**全 polynomial 恒等式**

```
F3(b²)(H,Y,W)=9216 b8(1+2b)² Fk(h,y,w)
```

を denominator を払って独立検証した。
k≠0,1,−1 は既 F3 の b≠0、1−4b²≠0 の unit chart に入り、
完全二次 Saito space は定数 GL2で一致、全六 quadratic magnetic pullback も可逆。
従って generic だけでなくこの unit chart の全 source は既 known103 内に入る。

k−1ではこの rational mapを使わない。
全 polynomial gradient matrix A(k) は generic rank16なので、全17minorは恒等零。
special k−1 で独立再生した16minor **−3913788948480≠0** は単元である。
従って近傍の全 kernel は正則 rank2 bundle。
free columnsを選び unit pivot inverseで basisを構成し、任意の special frameに定数 GL2で合わせて liftする。
Saito determinant はspecialで非零、よって近傍で単元。anchor及び commutatorも正則に解け、
任意の special磁気 R 全六係数をそのまま自由に保つ実 Poisson familyを作る。
known103の閉性で、**k−1の全 frame/全磁気も含む全 A3+A2 common-tangent source** が既成分に入る。
divisorの境界だけから異なる tensor の membership を逆向きに推論したものではない。

## prior111 と μ6九候補との全比較

MU6別査読で独立照合済みの immutable commit37d56e2a/blobdfc7c037 の全111行を再照合した。
同23 MPS29 は全 period actual SNC rank5 witness、ESSENTIAL5とGa9constantも actual essential5 witness、
新四 whole component は closed essential≤3なので相違。
MINIMAL_M の dense rank1 cofactor orbitに対し新全源はrank2。
旧 A0 は D5/H cusp tangent、新 D5はHが特異点を避ける smooth hyperflexであり、他三は異 profile。
旧 DPMは3A2/μ6、新全四はμ5。他の rank・dimension・closed multiplication/nonreduced・
旧 irreducible quintic witness の行も理由を保存し、全111件を処理した。

新四相互の全六組は quartic profile、同 A4+A1 の node/H incidenceで区別する。
μ6九候補との全36組は**同じ intrinsic quartic の totalμ が五と六で異なる**ため区別する。
MU6の36組と合わせ全十三候補の78組が揃う。
同23 componentなら dense23 actual orbit 同士が交わる必要があるが、intrinsic marked配置が異なり不可能。
この相互比較は μ6側の採録又は新component数を仮定しない。

root 四 ambient独立PASSと合わせ、四つの23成分採録を支持する。
A5 smooth-hyperflex と E6 smooth-flex の未知は保持し、他の E 残余をこの結果だけで閉じない。
