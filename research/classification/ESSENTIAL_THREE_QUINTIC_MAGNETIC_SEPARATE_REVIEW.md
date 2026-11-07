# 本質三変数既約quintic：全磁気除去・有限23軌道の独立review

2026-11-07。主義務E_reduced_nonnormal。**PASS**。
`ESSENTIAL_THREE_QUINTIC_MAGNETIC_REMOVAL.md` の全parameter・全Fの証明を別作者として査読した。
本質三変数の実既約被約quinticを持つ全rank四Poisson支持は、有限個のactual23次元PGL₅軌道の和へ帰着する。
具体的なcurve代表、既知成分への帰属、各閉包の全空間成分性は未解決であり、成分数は増やさない。

## 1. 実shearの符号と完全性

既受理のcanonical sourceを

`Π_R=U∧∂u+V∧∂v+[−(div V)u+(div U)v+R]∂u∧∂v`

とする。U,VはAR(F)の二次自由基底、Rは任意のS₂六係数である。
u'=u+l(a)、v'=v+m(a)は全parameterで可逆なunipotent linear変換である。
base括弧零によりanchorは変わらず、全normal括弧は

`R'=R+l div V−m div U−U(m)+V(l)=R+div(lV−mU)`。

これは近似でなく全tensor恒等式である。
独立checkerはU,Vの全三十六quadratic係数、R六、l,m六、GL₂四を自由にした**五十二parameter**で、
全十bivector coefficientのpushforwardとshearの両inverse compositionを展開して確認した。

## 2. kernel residue論証は全Fに適用する

δ(l,m)=div(lV−mU)=0としW=lV−mUと置く。
W(F)=0、div W=0、Wは次数三のbase場である。
θ=ιWΩは重み五の二形式なので、θ/Fは重み零、d(θ/F)=0。
ω=ιE(θ/F)はEuler-horizontal・重み零、closedなP²上rational一形式である。
そのnumeratorはΩ¹(5)のglobal sectionであるため、projective chartから余計な極は生じず、
全polesはFに沿うsimple polesだけである。

Fのgeneric smooth pointでω=(A dz+B dt)/zと書く。
closednessの二重極部分はBがzで割れることを強制し、単純極部分はA|Fの微分零を強制する。
Fは既約なので、このrational residueは一つの定数cである。
一般直線はFの特異点を避け、五smooth点で横断する。
そのP¹へrestrictしたωは他にpoleを持たず、全residueがcなので、residue定理により5c=0。
したがってc=0、Fのgeneric pointのpoleも消える。
他のdivisorial polesは無いので、smooth P²上の局所自由Ω¹のcodimension二延長からωはglobal holomorphic。
H⁰(P²,Ω¹)=0よりω=0。

ιEιWΩ=0はfraction field上EとWの平行性を意味する。
E(F)=5F≠0、W(F)=0からW=0。
AR自由基底U,Vの独立性からl=m=0。
従って六次元同士のδ:S₁²→S₂は**全admissible Fでisomorphism**である。

irreducibilityは不可欠である。可約Fではcomponent別の定数residueが次数付きで相殺し得るため、
この全六磁気除去を既済の可約quintic等へ無条件に転用しない。
logarithmic comparison theoremやFの特異点分類は、このkernel論証に不要である。

## 3. inverseとframe変換の全parameter復元

各Rに一意な(l,m)=−δ⁻¹Rがあり、actual shearがΠ_RをΠ₀へ送る。
磁気零は単なる軌道極限でなく全sourceの実同値normal formである。
Fを動かしてもAR₂はrank二kernel bundle。
各frame chartでδの六次正方行列は全閉点でdet非零なので、adjugate/determinant inverseはchart全域でregularとなる。
未知のspecial Fを除外したgeneric inverseではない。frame chartの補集合は他のframe chartで覆う。

vertexの実GL₂変換u'=Au+Bv、v'=Cu+Dvは

`U'=AU+BV, V'=CU+DV, R'=(AD−BC)R`

へ移す。この全十括弧恒等式とadjugate inverseも、上の五十二parameterのまま独立確認した。
AD−BC≠0はbasis frameの定義上の条件であり、degenerate frameはrank四の実Dを与えない。
よって固定projective Fの全frameと全磁気sourceは一つのactual PGL₅軌道に入る。

## 4. 有限curve軌道：constructible性の扱い

既受理の前稿で、全admissible Fは有限個のanalytic型multisetを持ち、
有限contact-determinacy jet-orbitのincidenceによるconstructible locus Xへ分かれ、各dim X≤8である。
ここでは**型multisetを固定したPGL₃不変のconstructible locus**を用いる。
任意のsmooth細分が群不変であると仮定する必要はない。

Saito自由基底の次数1,2,2からDer(−log F)₁=C E。
従ってprojective curve stabilizerのLie algebraは零、標数零なのでstabilizerはfiniteである。
全PGL₃軌道の次元は八。

各Xのclosureのirreducible componentsは有限で、次元は八以下。
八次元orbitのclosureはそのうち一つのcomponentに等しく、orbit自身はclosure内でopenである。
同じirreducible componentの二つのdisjoint nonempty open orbitは存在しない。
八未満の境界には八次元orbitは入らない。
従ってXのorbitは有限個、有限個のX全体でも有限個である。
これはdim≤8だけから有限性を推測したのではなく、全orbit dim=8とconstructible有限型を併用する議論である。
具体的なrational cuspidal quinticの分類リストは依存しない。

## 5. full PGL₅ stabilizerをcanonical liftで計算する

本質三変数空間は内在的に一意である。
実際、constant derivationでFを消す空間が二次元vertexで、そのannihilatorが本質変数Wである。
projective Dを保つPGL₅変換はWを保ち、base PGL₃に誘導するcurve stabilizerはfinite。
identity component内では全matrixのscalarを除き

`a'=a, (u',v')ᵀ=M(u,v)ᵀ+C a`

に限られる。このblock groupはGL₂四とshear六で次元十。

ここでprojective tensorの等号を勝手にraw liftの等号へ置換しない。
canonical divergence-zero liftはGL₅-equivariantである（constant volumeの変化は非零定数だけ）。
従って二つのprojective sectionがscalar倍で一致すれば、そのcanonical liftも同じscalar倍で一致する。
Euler ambiguityはdivergence-zero代表の一意性で消える。

Π₀のanchor U,Vの独立性から、これをprojectively保存するMはscalar tIだけである。
その後のshearが作るmagnetic項は、scaled frameに対応する同じ可逆δの像なので、C=0しかない。
逆にa'=a、u'=tu、v'=tvはΠ₀をt倍するため、全t∈C*がprojective stabilizerである。
identity componentは正確にC*、全stabilizerはそのfinite extension、dimension一。
actual orbit dimensionは24−1=23。

有限個のPGL₃ curve型と、各固定Fの全source一軌道性により、全Poisson支持は有限個の23軌道へ尽くされる。
その有限個の閉包は閉既約Poisson族であるが、より大きい成分の境界かもしれない。
次のE義務は有限代表の具体化とその全orbit closureの既知成分認識又は成分性である。

## 6. 固定証明書

finite job gate PASS後、52parameter全恒等式を0.027秒、最大RSS12,160KiBで確認した。
residue injectivity、有限orbit、stabilizerの証明は本文の数学的論証であり、sample計算で代用していない。

|artifact|SHA256|
|---|---|
|producer ESSENTIAL_THREE_QUINTIC_MAGNETIC_REMOVAL.md|bd255813e864fe6cfbb8a5f9c46b68e2bd5c9d4a78ae7c0218069bd34ff87a81|
|code/classification_essential_three_magnetic_review.py|7271da0dd5dc4b8a0bb255ed1e4592aea8709074874f4d7a79e8c5518ff0e35b|
|reports/CLASSIFICATION_ESSENTIAL_THREE_MAGNETIC_REVIEW_20261107.json|59f27fda831503ff8f2596248551ae708fa0425436fda3bf1f84397eb1fcad55|
|reports/CLASSIFICATION_ESSENTIAL_THREE_MAGNETIC_REVIEW_JOB_20261107.json|6788c7d10d9203e01e087c2e7bc4b38da134e5d72fbc13483ba2d7554fe4a90c|

全依存path/hashをjobとreportへ保存した。旧証明書、global台帳、原稿及びGitは変更していない。
