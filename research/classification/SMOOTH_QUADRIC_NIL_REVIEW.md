# 滑らかな二次曲面＋既約三次曲面：純冪零三型の独立除外

2026-11-07。課題 `E_reduced_nonnormal`。
標準発散零lift Πの実退化式をD=q c、qをrank五の二次式、cを既約三次式とする。
q-cofactor Xは非零なso(q)元であり、純冪零の可能性は(5),(3,1,1),(2,2,1)で尽くされる。
この入口は主ノート `SMOOTH_QUADRIC_CUBIC_SUPPORT.md` の既存証明を用いる。

## 1. 独立算術の射程

`code/classification_smooth_quadric_nil_review.py` は標準Fractionによる疎多項式演算だけを使い、
producer・SymPyをimportしない。
全150二次係数について、発散零、dq·Π=qNx、[N,Π]=0を全350係数slotに展開した。
表示sourceの全方向が厳密核に入り、独立であることと非零minorを組み合わせ、次のQQ階数を得た。

| 型 | 全線形条件の階数 | 全源次元 |
|---|---:|---:|
| (5) | 145 | 5 |
| (3,1,1) | 139 | 11 |
| (2,2,1) | 135 | 15 |

各affine source自身も全非斉次線形条件を満たすため、これで全線形源性が確定する。
計算は `jobs/smooth_quadric_nil_review.json` のgate後、約0.315秒でPASSした。
入力reportの全producer bytes/hash、独立checker、必要Jacobi係数のliteral位置、
minorを `reports/CLASSIFICATION_SMOOTH_QUADRIC_NIL_REVIEW_20261107.json` に保存した。
全十五変数のJacobi支持そのものは解いていない。除外に十分な条件だけを使う。

## 2. 型(5)と(3,1,1)

型(5)の全Jacobiには、非零定数倍を除いて次の四係数が実際に現れる。

    (n4−12)(n4−4),
    n2(3n4+4),
    16n2²+(9n4+12)n3,
    (15n4+20)n0+8n1n2+16n3².

従ってn4=4又は12であり、両方で3n4+4≠0なのでn2=n3=n0=0。
n1は任意で、両枝に代入すると全Jacobiが恒等零となることも独立確認した。
実退化式は順に `−x0 q²/2` と0である。
五つのPfaffian minorがこの式の勾配に一致することまで確認したので、
前者の三次因子は可約、後者はrank四ではない。

型(3,1,1)では、Jacobiを課す前の全十一源でEuler収縮Pfaffianが恒等的に

    x0 q Q2

の形になる。ここでQ2は座標次数二の式である。
Poisson点ではEuler収縮Pfaffianが実退化式になるため、
rank四ならc=x0Q2は可約である。全残Jacobiを解く必要はない。

## 3. 型(2,2,1)：Casimirからの必要条件

座標を `(u,a,v,b,w)=(x0,x1,x2,x3,x4)` と置くと

    N=u∂a+v∂b,       q=ub−av+w².

全十五源のbase(U=(u,v,w))上の三括弧は、厳密に次のJacobian bivectorである。

    {u,v}=F_w,  {u,w}=−F_v,  {v,w}=F_u,
    F=A(u,v)w+B3(u,v),
    A=(n8u²+2n12uv−n1v²)/2,
    B3=(n5u³+3n7u²v+3n10uv²−n0v³)/3.

F=0は七係数n0,n1,n5,n7,n8,n10,n12の消滅に等しい。
この代入で全五Pfaffian minorが恒等零になることを確認した。
従ってrank四ではF≠0である。

全35三次単項式に対するN作用は厳密に階数22であり、
その13次元核はbase三次式10個とq u,q v,q wで張られる。
よってN不変な三次式cは全て

    c=C3(U)+q ℓ(U)

と書ける。ここでℓは一次式である。
この核の完全性も別実装の非零minorと13個の独立核で確認した。

D=q cがΠのCasimirであり、NはUを固定するので、base関数との括弧から

    Q_base♯ dC3 + q Q_base♯ dℓ = 0

を得る。Q_baseの係数はUだけに依り、qはfiber変数a,bに非自明に依る。
多項式としてfiber係数を比較すれば、両項は別々に零となる。
例えばbの係数からu Q_base♯dℓ=0を得る。
これは整域での恒等式なので、u=0の点を除外する議論ではない。

Jacobian表示を使うと

    dF∧dC3=0,      dF∧dℓ=0.

FとC3は同じ次数三なので、第一式をEuler場で収縮すると
`F dC3−C3 dF=0`。F≠0より有理函数C3/Fの微分が零となり、C3=βF。
第二式でℓ≠0ならFはℓだけの函数で、同次性からF=αℓ³、α≠0。
この場合

    c=ℓ(αβℓ²+q)

は可約である。従ってcが既約ならℓ=0かつc=βFであり、Fも既約でなければならない。

## 4. 型(2,2,1)：必要Jacobi三係数だけでFの既約性を否定

全Jacobiのliteral係数には、非零定数倍を除いて次の三式が現れる。

    n1n8+n12²,
    n0n8+n1n7+2n10n12,
    n1n5−n10n8+2n12n7.

第一式はbinary quadratic Aの判別式である。従ってrank A≤1。
SL2は(u,v)と(a,b)へ同時に作用し、wを固定することでq,Nと体積形式を保つ。
この作用は全線形源を保つので、rank A=1ならA=αu²、α≠0へ正規化できる。
すなわち新しい正規化係数でn8≠0、n1=n12=0である。
残る二条件からn0=n10=0となり、Fはu²で割り切れる。
rank A=0ならF=B3(u,v)はbinary cubicであり、複素数体上で可約である。
いずれも前節の必要な既約性に矛盾する。

従って純冪零三型は全てD=q c、c既約の対象層から除外される。
正則半単純型の全MPS被覆と、階数二半単純・mixed型の独立除外を合わせると、
**滑らかな二次曲面＋既約三次曲面という因子型の全支持**はMPS40の有限閉和に含まれる。
この結論を他の因子型や全Poisson空間の尽くしへ拡張してはいない。
