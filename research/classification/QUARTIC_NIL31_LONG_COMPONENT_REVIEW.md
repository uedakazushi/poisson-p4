# nil31-long 第一五源：新22成分と既存96との差の独立査読

2026-11-07。`E_reduced_nonnormal`。
nil31-long 型5の第一完全五源を、別作者が全90係数・実PGL正規形・全周囲微分・
全96成分との相違まで独立に検証した。
**第一全五源の PGL 飽和閉包を、新しい22次元全空間成分として受理する。**
指定点は全 Jacobi scheme の滑らかな点で、この点を通る全成分は一意。
全五源は一つの実 PGL 軌道に入り、その軌道が成分内に稠密開である。
これは既存96に対する相対的な新成分性で、文献上の新規性は未確認。
全分類の上界や第二源の T40 帰属はこの査読の結論に含めない。

独立コードは `code/classification_quartic_nil31_long_component_review.py`、
報告は `reports/CLASSIFICATION_QUARTIC_NIL31_LONG_COMPONENT_REVIEW_20261107.json`、
有限契約は `research/classification/jobs/quartic_nil31_long_component_review.json`。
producer・producer helper・CAS・NumPy・legacy 実装の import はなく、
別作者が作った汎用 Fraction/sparse 算術だけを再利用した。
0.750秒、最大RSS20040 KiBで全検査 PASS。既存入力証明書 bytes は変更していない。

## 1. 全源と正当な適用範囲

元 cofactor は `Z=y∂x+z∂y+∂z`。
全90係数の div/Lie 全110行を独立に再構成し、階数75、全15自由係数を確認した。
全四 raw Jacobi、全係数表、実 Pfaffian と最高四次部分も独立演算と一致する。

全支持が二つの完全五源へ尽くされ、a=n14=0 の全補集合が H+既約quartic を持たないことは
既済独立証明 `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_REVIEW_20261107.json` に依存する。
この証明では第一源 n2=0、a非零を扱う。
自由係数は a,n3,n4,n8,n12 であり、残る全15置換を別に再入力し、元 source と一致させた。
これら全係数で strict Jacobi10式は零、h-Casimir斉次発散は

    N=y∂x+z∂y+h∂z

である。他の四自由係数の非零条件は一切ない。

## 2. 全五源を一つの実軌道へ送る正規形

canonical 値は a=1,n3=n4=n8=n12=0、すなわち n14=1,n13=−1/6、他n零。
その6倍整数代表は

    Pxy=3xz−2y², Pxz=3hx−yz,
    Pxw=12y²−18xz−wy, Pyz=2hy−z²,
    Pyw=6yz−18hx−wz, Pzw=6z²−12hy−hw.

canonical から全源への実 affine/PGL 変換は

    h'=h,
    x'=x−3n8 z−3a n4 w−(108a n4 n8+12n12 n4+2n3)h,
    y'=y−3n8 h,
    z'=z,
    w'=a w+(36a n8+6n12)h.

determinant は a非零。その逆も全係数で再構成した。
canonical tensor を実際に押し出し、全150係数が全源の二次斉次括弧と一致することを確認した。
単なる Pfaffian の一致や q 側だけの変換ではない。
tensor scalar factor は1。唯一反転する a の零補集合は既済の全支持証明で処理済み。
従って第一の全五源は、全てこの一つの実 projective PGL 軌道に入る。

## 3. 実quarticの既約性

unscaled canonical の実 Pf は F/6、上の整数代表では6Fである。ここで

    F=9h²x²−18hxyz+6xz³+8hy³−3y²z².

h=1 では多項式自動同型

    u=x−yz+z³/3, v=y−z²/2

により `F=9u²+8v³`。
u 二次式の判別式は C(v,z,w) の平方ではないので絶対既約。
F は h を因子に持たないので、その次数4斉次化も絶対既約である。
実退化五次式 χ は非零scalarを除いて hF。
従って被約で、唯一の一次因子 H=(h=0) と既約四次因子を持ち、一般Poisson階数は4。

参考として `F=−Disc(h,3z,6y,6x)/108` の全係数恒等式も再現した。
本質変数数を旧 nil4 との相違に使っていない。両者とも本質四変数である。

## 4. 全150変数・350式の独立成分証明

上の整数代表から、全350 projective Jacobi の全150方向微分 D を独立 Fraction 演算で作った。
25 Euler 垂直列、全25 GL5 列、scalar を含む G と、全五源の実微分を加えた Factual を作った。
全係数で DG=0、DFactual=0。独立な全有理消去・minor・全行列 hash を新報告に保存した。

    rankD=103, rank(vertical)=24,
    rankG=47, rankFactual=47.

従って全 C¹⁵⁰ の接次元は47、実族による局所次元下界も47。
局所環の次元と接次元が一致するので、その点は全 Jacobi scheme の regular point で、
通る全空間成分は一意である。
Euler24と射影scalar1を除くと成分次元は22、実 projective PGL 軌道次元も22。
この軌道の閉包が全空間成分そのものであり、全五源を覆う。
局所閉軌道の同次元性から稠密開軌道となり、一般PGL軌道次元22、Rosenlicht商次元0。
旧 nil4 の「接次元23、成分22」の場合と異なり、今回の滑らかさは全D103から直接従う。

## 5. 旧 nil4 の同じ1+4成分との内在的な相違

H が唯一の一次因子なので、H-Casimir homogeneous representative は内在的に一意である。
実際、二代表の差 E∧L が H-Casimir なら
`hLj−xjLh=0`。linear L について Lh=ch,Lj=cxj となり E∧L=0。
従ってその全斉次発散 N は PGL による共役と tensor の scalar 以外には依存しない。
H を定義する covector ℓ と N の dual-image flag の位置は内在的である。

今回の canonical は、並び(h,x,y,z,w)で

    Nnew=y∂x+z∂y+h∂z.

旧 regular-nil4 の既済22成分の実代表の発散は

    Nold=36(y∂x+z∂y+w∂z).

後者も旧点の全六括弧から新たに再計算した。
非零scalarを除くと双方 Jordan 型(4,1)、rank(N,N²,N³)=(3,2,1)、N⁴=0。
しかし ℓ=h について、全5×5有理行列の列空間を正確に計算すると

    new: h ∈ im(Nnew*)³,
    old: h ∉ im(Nold*).

第一全五源は上の実変換で canonical と同値なので、同じ flag を全域で持つ。
旧22成分も既済証明により22次元軌道閉包で、その canonical 軌道は稠密開。
もし両全成分が等しければ二つの稠密開軌道が交わり、この内在的 flag が矛盾する。
従って同じ1+4因子型でも異なる全成分である。
固定Jordan階数上の image flag を全空間の閉条件と仮定する論法は使わない。

## 6. 既存96成分の全監査

入力は `CLASSIFICATION_QUARTIC_NIL31_LONG_PRIOR96_20261107.json`。
全96 ID と次元・rank・一般因子型を固定した監査済み projection である。
既存各成分の成分性や一般因子型は、採用済み証拠を依存として使う。

* 一般rank2成分は閉rank条件で、次元22以外は今回の全成分性を証明した後の次元で区別する。
* 既存rank4・22次元成分は全15件。そのうち13件の一般五次式は非被約であり、
  全閉包が非被約五次式の閉条件に入る。今回の被約 hF の点はそこに入らない。
* `R4_TWO_HYPERPLANES_CUBIC_22` は全閉包が1+1+3の proper projective 積写像の閉像に入り、
  今回の既約四次因子を持つ点は入らない。
* 残る `R4_NIL4_HYPERPLANE_QUARTIC_22` は第5節の H–N flag で区別する。

非被約条件は `degG=1,2` の二写像 `(G,R)↦G²R` の閉像の和、1+1+3も proper な積写像の閉像。
退化式零の境界には affine cone に零を含めるので、全空間閉包での議論に欠落はない。
新報告に全96 IDごとの排除理由を記載し、未処理の同次元例外は残らなかった。

よって既存96に対して一成分を追加できる。第二源のT40帰属と合わせた型5全被覆は別途統合する。
