# AB は原点を含め全域で正規である

2026-11-02（Asia/Tokyo）。B=C[AB]、Dをその正規化とする。
既知のproper birational incidence f:X→AB と、
Z=LG(2,4)へのaffine射影π:X→Zを用いる。

## 1. 結論

    B=D、すなわち被約ABは原点を含め全域で正規である。
    Rⁱ f_* O_X=0  (i>0).

従来は原点以外の正規性と次数三までのB_n=D_nが既知だった。
新しい四次の等号だけから全次数を外挿するのではない。
以下でDのS-module生成次数が4以下であることを別に証明し、有限個の等号からB=Dを導く。
ここでS=Sym(W45*)であり、各二次テンソル係数の次数を1とする。

この結論は被約ABについての主張である。実際のT01の零一次ファイバーの厚み、
全S61の成分の尽くし、成分数の上界は解決していない。

## 2. incidence の局所graded模型と低次Tor

`AB_INCIDENCE_20261101` により、XはZariski局所で

    U × Y × A³_g,    Y={N∈sl3:rank N≤1}

に同型で、Nとgに次数1を与えるgraded同型である。
補Lagrangianの変更はNを保ちgへNの線形式を加えるので、gradingを保つ。
Nの八座標とgの三座標からQへの線形写像はrank11である。
既知の45個の独立テンソル係数の制限でもrank11であることを再検算する。
従ってS⊗O_Uの45一次生成元からこの11一次座標への写像は全射で、
そのkernelはrank34の局所自由層である。局所的に線形に分裂させられる。

Yの八変数多項式環上の最小自由分解は、rankと内部次数が

| homological degree k | rank | generator degree |
|---:|---:|---:|
| 0 | 1 | 0 |
| 1 | 9 | 2 |
| 2 | 16 | 3 |
| 3 | 9 | 4 |
| 4 | 1 | 6 |

である。一般traceless 3×3行列の全2×2 minorsからQQ上で再計算し、
全正homologyの消滅、最小性、全shiftを検査する。
三つの自由g変数を追加してもこのTorは変わらない。
さらに34個の零となる一次変数のKoszul複体をtensorすると、
Sに関するTorの次数は局所的にこの表と34一次Koszulの畳み込みになる。

m=S_+、k=S/mとし、

    T_i = Tor_i^S(O_X,k)

と置く。これらはmに殺され、scheme論的零ファイバーZ上の局所自由層である。
その内部次数は、i=0なら0、1≤i≤3なら高々i+1である。
具体的にはi=1,2,3でそれぞれ

    {1:34, 2:9},  {2:561, 3:322},  {3:5984, 4:5602}

となる。大域的なparabolic extensionの分裂や、大域的な相対最小分解の存在は仮定しない。
局所Torの次数消滅だけを使う。

## 3. 全次数でH²以上が消える

π_*O_X=⊕_(n≥0) E_nと書く。
既存二次・三次・四次ノートの、垂直g次数jによるP-stable filtrationを任意nへ適用する。
associated gradedの各部分はY_(n−j)⊗Sym^j(V_vertical)である。
coreの一次Levi weightsの和は2、verticalの和は4なので、
そこに現れる全てのLevi最高weight(a,b)は

    a≥b,   a+b=2n+2j≥0

を満たす。Borel–Weil–Bottでρ=(2,1)を足しu=a+2、v=b+1とすると、

    u>0,   u−v=a−b+1>0,   u+v=a+b+3>0.

従ってu>|v|であり、vが正なら長さ0、負なら一回の符号反転で長さ1、
v=0ならacyclicである。次数2以上のcohomologyは生じない。
有限filtrationを戻しても

    H^q(Z,E_n)=0  (q≥2, n≥0)

となる。これは有限個のnの表に基づく推測ではない。
E_n=0(n<0)も用いる。πがaffineで、Zの有限affine被覆でcohomologyとgraded直和が可換なので、
K=RΓ(X,O_X)はH^q(K)=0(q≥2)を満たす。

## 4. 導来ファイバーから生成次数を制限する

fをSpec Sへの射とみなす。proper性からKの各cohomologyは有限graded S-moduleであり、
cohomological degreeは0から3に限られる。
Sは多項式環なので、Kは有界な有限graded自由複体で表せる。
同次数の単位成分を順に消去し、全微分の成分がmに属する最小複体Fを取る。
すると

    F^q⊗_S k = H^q(K⊗^L_S k)

である。最小化はgradingを保つ。

kはS-perfectなので導来projection formulaを適用できる。
その右辺をTor層から計算するhypercohomology spectral sequenceは

    E2^{j,−i}=H^j(Z,T_i)  ⇒  H^{j−i}(K⊗^L_S k)

となる。dim Z=3ゆえj≤3である。
total cohomological degree1にはi≤2だけが寄与し、§2から内部次数は高々3である。
total degree0にはi≤3だけが寄与し、内部次数は高々4である。
全differentialは内部次数を保つので、subquotientへ進んでもこの上界は変わらない。
従って

    F¹の生成次数≤3、F⁰の生成次数≤4.

同時にF^q=0(q>3)である。ここではTor spectral sequenceを退化させてはいない。
大域的にTor束を直和分裂させる必要もない。

## 5. 高次cohomologyを消し、全次数生成性へ進む

§3によりH^q(F)=0(q≥2)である。
もしFに次数2以上の非零項があれば、その最大次数rを取る。
H^r(F)=0ならF^(r−1)→F^rは全射であるが、最小性からmで割った写像は零であり、
非零自由加群F^rへの全射にはなれない。graded Nakayamaに矛盾する。
従ってF^q=0(q≥2)である。

ゆえにH¹(K)=coker(F⁰→F¹)は次数3以下で生成される。
既存証明からH¹(Z,E_n)=0(n=0,1,2,3)である。
n=0はO_Z、n=1は(3,−1)のacyclic性と(2,0),(3,1)の高次消滅、
n=2,3はparabolic評価を用いた既証明の消滅である。
負次数も零なので、この有限生成上界によりH¹(K)=0が全次数で従う。
同じ最小性とNakayamaによりF¹=0になる。

従って

    D=H⁰(K)=coker(F^(−1)→F⁰)

はS-moduleとして次数4以下で生成される。
新しい四次証明を含め、B_n=D_n(0≤n≤4)であるため、
全生成元をB内に選べる。BはSの像なのでS·B⊂Bであり、D⊂Bとなる。
逆包含は既知だからB=D。Xの正規性とproper birational性からDはBの正規化であり、
ABは原点を含め正規である。
同時に全H^q(K)(q>0)が消え、affine底上のR^qf_*O_Xも零となる。

## 6. 基礎文献・再現・限界

導来議論の基礎は、2026-11-02に本文を確認した次の一次資料に従う。

- [Stacks Tag 02O3](https://stacks.math.columbia.edu/tag/02O3)：proper射のcoherent高次像。
- [Stacks Tag 066Z](https://stacks.math.columbia.edu/tag/066Z)：regular環上の有界有限cohomology複体のperfect性。
- [Stacks Tag 00MR](https://stacks.math.columbia.edu/tag/00MR)：自由複体の単位成分の消去。
- [Stacks Tag 01E6](https://stacks.math.columbia.edu/tag/01E6)：perfect因子に対する導来projection formula。

上の生成次数基準とABへの適用は、これらの基礎を組み合わせた本ノートの議論である。
局所Torのshift、rank11、次数0..4の既存証明との整合を専用producerで確認する。
独立概念レビューは別ノートに記録し、プログラムだけで導来議論が証明されたとはしない。

実際のT01零ファイバーの非被約構造、全S61の尽くし、一般アフィン五点以外の不存在、
全成分交叉・軌道分類は別の未解決課題である。
