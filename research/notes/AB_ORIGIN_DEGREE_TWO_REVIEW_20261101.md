# AB 正規化環の次数二生成性の独立点検

2026年11月1日。[主証明](AB_ORIGIN_DEGREE_TWO_20261101.md)の幾何学的議論を別担当者が点検した。
結論 D2=B2、dim D2=761 を受け入れる。原点の正規性や全次数生成性は含まない。

## 点検した論理

proper birational な正規 incidence X→AB と AB の affine 性から、正規化環は H0(X,O_X) である。
各次数の fiber coordinate ring は一定次元の P-加群であり、LG 上の同変局所自由層になる。
core の座標は下三角 symplectic shear で変わらず、vertical の座標には core の一次式が加わる。
従って vertical 因子の個数による filtration は安定であるが、一般には分裂しない。

主証明の次数二の Levi 分解と9本の独立二次 minor の扱いをコードで確認した。
Borel–Weil–Bott には [Fonarev, Theorem 3.2](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9880&wshow=paper)
の一次資料本文を読んだ。Σ^(a,b) U* には同定 Σ^(a,b) U*=Σ^(−b,−a) U を用い、
定理の −lambda+rho が (a+2,b+1) となる。従って各 graded piece の (h0,h1) は
(84,84)、(486,0)、(275,0) であり、唯一の高次コホモロジーは Sym6(V) の84次元である。

接続写像が非零であることについて、次の二点を区別して確認した。

1. 非零 G-部分加群の固定平面での評価は零写像ではない。もし全部零なら、
   G の推移性と部分加群の G-安定性により、全ての点で全ての切断が零になるためである。
   従って Sym6(V) が H0(E2) に現れれば、非零 P-写像が fiber へ得られる。
2. y の P 軌道は p 成分が非零な一次形式の開集合を含む。GL2 で任意の非零 p 成分を取り、
   symmetric shear で任意の q 成分を加えられる。よって y6 は P-cyclic である。
   これは「G-既約だから P-cyclic」と推測したものではない。

独立コードは parabolic の7個の行列から Sym6 の作用を再構成し、y6 の巡回空間の階数84を厳密に確認した。
fiber の weight (−6,0) は core×vertical の h0g0 だけである。
core² と vertical² はそれぞれ中心 weight の合計が −4、−8 で、この weight を持たない。
主証明と同じ pullback convention を両方の加群に一貫して用いればよく、逆作用に揃えても議論は変わらない。

さらに主コードの tensor extraction を使わず、q,p 順の block 行列を直接掛け、
w→w+u z が h0 を固定し g0→g0+12u h0 とすることを確認した。
この shear は y6 を固定するが h0g0 の微分は12h0²である。
h0² は core の rank-one traceless cone 上で零でない（N=E13）。
従って equivariant な像 c h0g0 は c=0 を強制され、cyclicity に反する。
このため Sym6(V) は H0(E2) に現れない。

core×vertical の H0 と core² の H1 の間で可能な接続写像は、唯一の Sym6(V) 同士の写像である。
もし零なら半単純性によって Sym6(V) が H0(F1) に持ち上がり、H0(E2) にも現れてしまう。
従って接続写像は84次元の同型である。長完全列より h0(E2)=845−84=761、H1(E2)=0 を得る。

## 係数計算の再現と範囲

主スクリプトを新しい出力先へ再実行し、有限体32003上の選択761次 minor の積11854を再現した。
コードの denominator 検査・疎行消去・pivot 保存を読んだ。
この有限体計算の役割は QQ 上の下界761だけであり、上界は上の特性零の幾何学的議論から得る。
二つが一致するので、実際の45個の次数一係数の積が正規化環の次数二全体を張る。

別コード `code/ab_origin_degree_two_review_20261101.py` と
`reports/AB_ORIGIN_DEGREE_TWO_REVIEW_20261101.json` は、巡回性・block shear・weight・Bott 次元を独立に確認する。
761次係数行列を別実装で再構成したとは主張しない。

```
.venv/bin/python code/ab_origin_degree_two_review_20261101.py --output /tmp/ab-origin-degree-two-review.json
```

原点外の正規性や次数二までの一致から全次数生成性は従わない。
残る問題は次数三以上の追加正規化関数の有無である。
