# 三重直線 m二 constant-double：全八行列と全P2を保つ二nil源の独立確認

2026-11-08、PRIMARY G/E、別作者 `/root/quartic_review`。
**六 semisimple 全源は空、二 nil 全源は空又は既 essential-three 八成分へ全点被覆される**ことを独立受理する。
対象は先に独立確認済みの R=3L、normal 位数二、constant double root の全八 first symbols。
他 root degrees、SF、R=2L+M、m三以上を含めた全分類は未結論。

## 1. 二つの旧欠陥を保持して訂正する

旧 producer report `CLASSIFICATION_PLANE_TRIPLE_LINE_M2_DOUBLE_EIGHT_FULLHIGHER_20261108.json`
(`5dbe47a0b0ec2266c7809bf36f817e057fdd6b9b20a52685ecebf0fb08cdf161`) は不変である。
別作者として次の二点を指摘し、root は明示訂正 archive と新 V2 report を保存した。

* 全行の `actual_P0_zero=true` は誤 metadata。実計算は非零 `actualG2` と一致しており、その方を検証する。
* NIL0 の旧 a-factor 判定は homogeneous K2 basis のみを用いたため無効。
  固定 P1 の P2 は J2 P2=−Jac(P1,P1) の affine torsor の部分集合であり、一般に K2 の元ではない。
  P0+K1+K2 の a-row だけから、全解の a-row を a 倍とは言えない。

この二点を受理根拠から除いた。旧 NIL0 全空は撤回され、以下の Ore 証明による全点被覆へ変わる。
V1/V2 の全150P0、350×30/9、150×39 group、全 kernel/minor、coarse bound は全て同一と独立確認した。
V2 の修正は上記 metadata と full9 P2 へ戻した factor flags である。
旧 source/report/job は correction archive に exactbytes が保持されている。

## 2. 全八行列の別 Fraction 実装

既独立 checker へ isolated `--triple-line-double-eight` だけを追記した。
旧 concurrent17 code `605cfaa0...`、当時 job/manifest の全文は公開 preappend snapshot
`plane_m2_concurrent_reviewer_605cfaa0_preappend_20261108.json` に保存し、旧17 FunctionDefs の AST は全て不変。
旧成功 report は変更せず、current old job は明示 compatible-extension 理由を持つ。
producer/CAS を import していない。

有限 gate 後 **0.816秒 / 22,328KiB**、独立 sparse/Fraction run は全 PASS。
全150P0 と非零 actualG2、全350 Jacobi、350×30 J1、350×9 J2、150×39 Euler/Levi/scale、
complete kernels・全行代入・非零 minors を元の定義から再構成した。
K1 は unconstrained 全核であり、m三の leading-target 条件を誤適用していない。

| profile | dim K1 | dim K2 | firstsymbol affine B | coarse projective upper |
|---|---:|---:|---:|---:|
| sem ±1/2, ±1 | 9 | 0 | 9 | 23 |
| sem ±2 | 8 | 1 | 9 | 23 |
| NIL0 | 12 | 3 | 8 | 28 |
| NIL1 | 11 | 3 | 9 | 28 |

六 sem は全 K1 bad u-row 零である。**全九 P2** に normal column が無いので全 actualD は u 因子を持ち、指定既約源は空。
両 nil の full9 での u/a factor は共に false、粗上界28はそのまま未閉として保持した。
以下はこの粗上界を小さい値へ書き換えたものではなく、全式からの別構造証明である。

## 3. 全 K1 と全九 P2 の Euler regauge を独立再構成

producer 新純証明 `PLANE_TRIPLE_LINE_M2_CONSTANT_DOUBLE_ORE_WHOLE_BOUND.md`
(`ed1ae0afb6f07e4b3daaf7c50f6b423e457198fec6e58788f0390ebb78b7190e`) を独立に全読解・再導出した。
両 nil の full K1 basis の c 依存は Pbc の c*u/c*v だけ。
B=Pbc_c とし、T=B∂b、P'=P+E∧T とする。
Pbc へ −cB、Pab へ aB、Pbu へ −uB、Pbv へ −vB が加わり、全 c 依存が消える。
これは同じ projective tensor の Euler gauge であり、有理座標同値ではない。
全 P2 九は元から c 非依存で、ab の全 arbitrary F2(u,v) は Qab、ac/bc の六係数は W に入る。

独立 checker は各 K1 basis と各 **full P2 basis** にこの操作を行い、producer の NIL0 十二係数/
NIL1 十一係数の Q、B の全線形式と一つずつ一致させた。
従って各 homogeneous kernel のみならず、任意 affine particular を含む全 P2 を保持して

    P'=Q(a,b,u,v)+∂c∧W(a,b,u,v)

が必要である。W を固定せず、元 source との全相関は必要包含として残す。

projective Jacobi の c 係数は ∂c∧[Q,Q] なので strict [Q,Q]=0。
残定数係数は E_y∧[Q,W]=0、y=(a,b,u,v) を与える。
元 P の strict Jacobi を仮定していない。

## 4. constant factor と actual quintic の全点認識

Q=∂b∧V、V が b 非依存なら PfQ=0。
従って actualF は c 非依存、sign を除き det(E_y,∂b,V,W) である。
module E_y∧[Q,W]=0 に ∂b を wedge すると
det(E_y,∂b,V,∂bW)=0。
det の b 微分では ∂bE_y=∂b の項は重複で零、∂bV=0 なので、これは ∂bF=0 に等しい。
ゆえに actualF は b,c の両方に依存せず essential 高々三。
非零既約 homogeneous quintic は essential 二以下ではあり得ないため、既 accepted
essential-three 全八成分 theorem を適用できる。
W の source dimension、homogeneous/particular の区別を追加仮定しない全点命題である。

NIL0 の独立抽出 Q は常に ∂b∧V、b 非依存である。
従って NIL0 は全空とは言わず、全 rank四・既約 actualD 解をこの八成分へ認識する。

## 5. NIL1 の二 strict Jacobi 係数

独立抽出 Q は producer 記法で

    Qab=a((2A+3T+J)u+(10m/3)v)+F2(u,v),
    Qav=T*u²,
    Qbu=−A*u²−(14m/3)u*v,
    Qbv=a*u+R*u²+(J−A)u*v−(23m/3)v²,
    Qau=Quv=0.

元の strict cyclic definition から手で直接再計算すると

    Jabu=−(14/3)T*m*u³,
    Jabv=−T(A+3T)u³−(28/3)T*m*u²*v.

F2 の全三係数は Quv=Qau=0 と b 非依存により両式へ寄与しない。
従って全複素解で Tm=0、T(A+3T)=0。
T=0 では Q が再び ∂b∧V なので前節の全 essential-three 認識を適用する。
T!=0 では m=0。
元 K1 の唯一 bad u-row は m*v²*∂c であることも各 basis から独立に確認した。
Euler gauge は u-row に u 倍しか足さず、full9 P2 に normal column は無い。
よって m=0 の全 full solution は u 因子を持ち、actualD 既約性と矛盾する。
T/m/A の零境界、F2 の全係数、affine particular を落としていない。

これで二 nil 全源も空又は essential-three 八成分へ入り、先の全八必要 first-symbol cover と合わせ、
**R=3L・m二 constant-double 全源の pointwise E 被覆を受理する。**
他 flags/R/m 及び全分類は別残余である。
