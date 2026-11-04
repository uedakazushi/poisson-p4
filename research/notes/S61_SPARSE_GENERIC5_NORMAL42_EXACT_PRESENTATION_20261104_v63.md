# 実五変数源の閉残余を42法線の次数四以下idealで表す

2026-11-04。既存 `S61_AFFINE_NORMAL_FINITE5_20261103_v63` に同型のscheme elimination定理がある。本稿はその抽象型の新規性を主張せず、今回独立確認された実K5源・130残差へ、固定126minorを仮定しない**全入力付きtyped operator**を与える。大きいGram係数や逆行列は計算・展開しない。残余の空性は未証明である。

## 1. 元のschemeとcanonical parameter左逆

K=QQ(H,L,s0,q0,t0)、R=K[n0,...,n41]とする。選択五式をF_r(z,n)、元129残差と実G6残差をR_i(i=1,...,130)と書く。これらは元134式と同じidealを定め、G6を含めても新しいアフィン条件を加えない。全式の未知変数次数は2以下で、全法線二次項を保持する。

u_hk=z_h*n_kを126列に並べる。受理済みの実回路から

```
R130=X130x126*u+b130(n)
```

を取る。b_iは実R_iのz=0への特殊化で、次数≤2、b_i(0)=0である。Xはparameterだけの行列でrank126と証明済み。

Kはformally realな体である。Cauchy–Binetにより
det(XᵗX)=Σ_I det(X_I)²≠0
だから、Lを唯一のparameter行列

```
(XᵗX)*L=Xᵗ,  L126x130
```

として定められる。LX=I126である。全係数を固定した有限線形系であり、未知変数を反転する手続きではない。typed recipeの行列積は有限和として定義するが、126×126 Gramの2063880個のscalar積を列挙しない。Gram inverseや126 Gaussianを実行したとの主張はしない。

𝒫(n)=−Lb(n)、H(n)=b(n)+X𝒫(n)と置く。𝒫は126個の次数≤2のnormal多項式で、𝒫(0)=0。H130も次数≤2で、係数行空間は証明済みの四次元pure-normal空間を張る。四本の基底を抽出せず130個をそのまま用いる。

形式的実性を使うのはK上でGramが非零なことだけである。Cへの係数拡大後も非零detは非零で、複素解nに順序や正値性を仮定しない。Σn_j²を反転する操作は一切用いない。

## 2. 2923個のnormal方程式

J⊂Rを以下で生成する。

1. H_i(n)の130式。次数≤2、線形包は四次元。
2. 全h=0,1,2、j<kについてC_hjk=n_j𝒫_hk−n_k𝒫_hj。2583式、次数≤3。
3. 全選択五式r、全jについてQ_rj=n_j²F_r(𝒫_0j/n_j,𝒫_1j/n_j,𝒫_2j/n_j,n)。210式、次数≤4。

最後の表示は**定義時にn_jを反転しない**。F_rのplane、linear-normal、cross、normal² blockを使い

```
Q_rj=n_j²*F_r(0,n)
     +n_j*Σ_h 𝒫_hj*(coefficient_of_z_h_in_F_r)
     +Σ_(h≤l) coefficient_of_(z_h z_l)*𝒫_hj*𝒫_lj
```

として多項式を組み立てる。F_rの全次数が2以下、𝒫の次数が2以下なので次数は4以下になる。全45未知変数における次数≤2の1081単項式を全42chartで扱い、このhomogenization恒等式45402件を確認する。normal²項はF_r(0,n)に含め、捨てない。

したがって、ここで定めるJは130+2583+210=2923生成の冗長な実idealである。四本のpure-normal基底を選べば2797生成にできるが、その基底を計算したとは主張しない。

## 3. 各D(n_j)でのscheme同型

A=K[z,n]/(F5,R130)とする。A内でR=0とLX=Iからu=𝒫、H=0を得る。従って全CとQも零で、R/JからAへの自然写像がある。

逆に(R/J)_(n_j)上でz_h=𝒫_hj/n_jと定める。Cの等式から、すべてのkについて𝒫_hk=z_h n_kとなる。よって

```
R=Xu+b=X(u−𝒫)+H=0.
```

またQ_rj=0とn_j²の可逆性からF_r=0である。この二つの写像はnでは恒等であり、zもu=𝒫から戻るので互いに逆。従って

```
A_(n_j) ≅ (R/J)_(n_j)
```

はidealと商環についての厳密な同型である。radicalだけの同型ではなく、全nilpotentと埋め込み構造を保つ。chart間でも同じcompatibilityがzの一致を保証する。42個すべてのD(n_j)を使い、一つのΣn² chartで代用しない。

## 4. 有限性、閉埋込み、原点を含めない降下

cross rewrite z_h n_k=𝒫_hk(n)はnを含む単項式のz次数を減らす。既知plane五二次式は最高二次rank5、その一次倍は全三次rank10である。G6 packetのcentered triangle証明からこのことは今回のK上でも成立する。従って純z三次も次数≤2へ戻り、R上のAの生成数は高々五である。これはR有限性でありK有限性ではない。

既存の実polarization splitterによってA=E5×B、かつA/(n)=E5である。従ってBはR有限、m=(n)でB/mB=0となり、中山の補題でB_m=0を得る。よってfiniteなnormal射の像は閉で、原点を含まない。

一方D(n_j)上ではE5が消え、前節からR_(n_j)→B_(n_j)は全射である。R→Bのcokernelは有限R加群で、mではB_m=0、他のすべてのprimeではいずれかのn_jが単元なので零になる。従ってcokernelは零で、R→Bは全射、normal射は**閉埋込み**である。

この有限cokernelの議論が、open chartの表示を実際の閉残余全体へ降ろす要点である。chartの点集合が一致するだけで閉埋込みと呼んでいない。

## 5. 正確なsaturation idealと空性criterion

I_B=ker(R→B)とする。各D(n_j)でI_BとJは等しく、Spec Bは閉で原点を避けるので

```
I_B = intersection_(j=0)^41 (J:n_j^infinity)
    = J:m^infinity,
B ≅ R/(J:m^infinity).
```

これはschemeの等式である。有限の通常多項式eliminationとしても

```
J:m^infinity = intersection_j ((J,1-t_j*n_j) intersect R)
```

と定義でき、42chartをすべて保持する。したがって

```
B=0 ⇔ J:m^infinity=R ⇔ あるdでm^d⊂J.
```

今回このcolon、冪所属、unitの係数計算は行わない。typed recipeが定める全2923式とその入力は保存し、次の尽くし証明が検証すべき有限idealを固定する。

## 6. 保存と旧成果との差

同じscheme-elimination型は旧NORMAL_FINITE5にあり、当時は特定126minorと129残差、三本のpure-normal関係を用いた。今回の保存対象は、実五変数源に束縛した130残差・四次元関係空間・canonical Gram operatorと、全2923生成のtyped polynomial constructionである。old107MB cross係数や固定旧minorを今回の源で非零だと仮定しない。

`code/s61_sparse_generic5_normal42_exact_presentation_20261104_v63.py --output work/NEW-output` で一回60秒・2GB・32MiB。明示した新parameter多項式supportは0、Gramはlazy matrix operatorのまま。新しい点、素数、GB、Gaussian、fullcomponent尽くしの計算は使わない。

## 確定した構成結果

一回の構成・検査は0.7035秒で完了した。保存operatorは24501 node、literal JSONは3164022 byte、gzip約144KiBである。新しい明示parameter多項式項数は0。Gramを展開すれば必要な2063880 scalar積は一つも列挙しておらず、parameter線形系の係数解や逆行列も評価していない。

全2923生成がK[n42]の多項式として型付けされ、130個のHは次数≤2、2583個のCは次数≤3、210個のQは次数≤4と確認した。すべてnormal原点で零となる。全42×1081=45402個のhomogenization単項式確認も完了した。parameter solveの入力はMatK型の行列だけで、nやzを含む型の逆元は構成に存在しない。

`typed-normal42-operator.json.gz` は全source root、X、b、Gram、Lの定義、𝒫、H、全compatibilityとquartic、全42chart写像、saturationの有限elimination recipeを含む。`scheme-elimination-theorem.json` は実A、exact B、閉埋込みとunit criterionを区別し、saturation未実行を明記する。既知五点clopen分離については元のpolarization primaryと独立reviewを入力hashへ追加して固定した。この最後のhash/gate追記は依存の明示だけで、算術構成を再実行していない。

したがって今回確定したのは、次数四以下の有限ideal Jと**全schemeを保つexact presentation**である。Jのscalar係数全展開、J:m∞の生成元計算、mの冪の所属係数、残余の空性はいずれも次の仕事として残る。
