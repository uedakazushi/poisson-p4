# rank三二次因子と三次因子：Anil3 の全五源

2026-11-07。課題 `E_reduced_nonnormal`。
`DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTION.md` §5 の五型
`(Anil3,Cnil2)` の upper/lower/zero と `(Anil3,C0)` の nonzero/zero を扱う。
全150係数の線形源と全Jacobi係数は
`CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json` の独立再構成を入力とする。
全五源の既約三次因子を持つ被約支持は **三つの一軌道族**へ帰着する。
§6の忠実graph認識により upper はT39に、§7の全パラメータ実曲線により
残るI・IIも既知MPS成分の和に含まれる。従って対象の全五型の被覆が完了する。
ここで新しい成分は加えない。

記号は `(a,b,c,r,s)=(x0,x1,x2,x3,x4)`, `q=2ac+b²` とする。
以下の `C` は Euler 収縮 Pfaffian を q で割った多項式である。
Jacobi を課した後だけ、これを実三次因子と同一視する。

## 1. Cnil2 upper：全二源だけが残る

cofactor は `N=a∂b−b∂c+c∂r+r∂s`。
全七変数源の必要Jacobi係数には

    (n6−2)(8n6−11), n5(n6−1),
    5n2−8n4n6+16n4,
    (24n6−33)n1+(6n6+48)n3+8n5²,
    (−9n6+8)n0−2n1n5+3n2²+9n2n4+4n3n5+6n4²

がある。従って n5=0 で、全支持は次の二枝で尽くされる。

- `n6=11/8, n2=−n4, n3=n0=0`。全Jacobiは零だが C=0。
- `n6=2, n2=n5=0, n1=−4n3, n0=3n4²/5`。n3,n4 は任意。

第二枝の三次式は

    C=C0−n3*a*q/2−n4*(a²r−abc−b³/3)/2−n4²*a³/10,
    C0=5(12acs−9ar²+6b²s+6bcr+2c³)/72.

`∂C/∂s=5q/12`、`C|a=b=0=5c³/36` である。
q は既約なので、s に関する一次多項式 C の係数と定数項は互いに素。
Gaussの補題により、全二源で C は既約である。

基点 n3=n4=0 から任意の源への射影線形変更は

    a'=a, b'=b, c'=c,
    r'=r−2n4*a/5,
    s'=s−2n4*b/5+6n3*a/5.                         (1)

行列式は1で、q,Nを保つ。
基点の全十括弧を (1) で pushforward すると任意の源の括弧になることを、
全二自由係数のまま厳密に検証した。従って全二源は一つの PGL 軌道である。

## 2. Cnil2 lower と zero は空

### lower

cofactor は `a∂b−b∂c+(r+c)∂s`。
必要Jacobi係数

    (2n2−1)(4n8+3), (n8+2)(n2+3n8+3), n7(2n2+6n8+1)

から `(n2,n8)=(1/2,−2),(1/2,−7/6),(−3/4,−3/4)` かつ n7=0。
残る係数を解くと、それぞれ

    (n4,n5,n3)=(-4n6,0,27n6²/5),
    (n6,n4,n3)=(0,3n5,9n5²/5),
    (n6,n4,n3)=(0,n5,0).

いずれも他の全自由係数を残したまま、全Jacobiは恒等的に零となる。
第一枝の C は

    −(6n6*a+5r)(18n6*a²+10q+15ar)/120.

第二枝では a が C を割り、第三枝では C=0。
従って既約三次因子を持つ全支持は空である。

### zero

cofactor は `a∂b−b∂c+r∂s`。
全十一源のJacobi係数に 5n10/4 と n5² があるので、被約支持で n10=n5=0。
この代入後、全 C は a で割り切れる。他の九変数の制約を解く必要なく、対象支持は空。

## 3. C0 zero は空

cofactor は `a∂b−b∂c`。
全十二源のJacobi係数に n10*n9 と n11*n9 がある。
元の全 C に対し、a=0 で

    C|a=0=b²(n10*n9*r+n11*n9*s)/5.

従って全Poisson支持では a が C を割る。対象支持は空である。

## 4. C0 nonzero：二つの全四源

cofactor は `N=a∂b−b∂c+c∂r`。
全九源の必要Jacobi係数は、非零定数倍を除いて

    (n8−2)(2n2−2n8+3), n7(3n2−3n8+2),
    2n1n7+(6n2−3n8+3)n5,
    (2n2−2n8+3)n4+(2n2−6)n6,
    6n0n7+(18n2−12n8+10)n3+3(n4−2n6)n5.        (2)

### I：n7≠0

最初の二式は n8=2,n2=4/3 を強制する。
残る三式より

    n4=2n6, n5=−2n1n7/5, n3=−3n0n7/5.             (3)

自由係数は n0,n1,n6,n7、開条件は n7≠0。
全元Jacobiはこの全四源で恒等的に零となる。

    C=(6n0n7*a³+6n1n7*a²s−30n6*a*q
       +30n7*a²r−30n7*abc−10n7*b³−25q*s)/60.

r に関する係数は n7*a²/2、a=0 の定数項は
`−b²(2n7*b+5s)/12` なので、全開条件で C は既約である。
基点 `(n0,n1,n6,n7)=(0,0,0,1)` から全源への変更は

    r'=r+A a+B s, s'=δs+C a,
    δ=n7, C=−6n6/5, B=−n1n7/5,
    A=−n0/5+6n1n6/25.                              (4)

a,b,c は固定する。行列式は δ≠0。
全十括弧の恒等式を検証したので、全源は一つの軌道である。

### n7=0 の補集合を保持する

n8≠2 なら (2) から n2=n8−3/2 であり、この代入だけで a|C。
n8=2 では残る必要式が

    (2n2−1)n5=0,
    (2n2−1)n4+(2n2−6)n6=0,
    (18n2−14)n3+3(n4−2n6)n5=0

となる。
n2=1/2 なら n6=0,n3=3n4n5/5 で
`C=−n5*a²(2n4*a+5s)/20`。
n2≠1/2 なら n5=0。さらに n3=0 の全枝では q|C。
従って既約 C を持つ残余は n3≠0 の場合だけである。

### II：n7=0,n3≠0

上の式から

    n8=2,n2=7/9,n5=n7=0,n4=8n6.                  (5)

自由係数は n0,n1,n3,n6、開条件は n3≠0。
全元Jacobiは恒等的に零であり

    C=−[6n3*a³+(18n6*a+5s)q]/36.

s に関する係数 q と定数項は互いに素なので、C は全開条件で既約。
基点 `(n0,n1,n3,n6)=(0,0,1,0)` から全源への変更は (4) と同じ形で

    δ=n3, C=−18n6/5, B=−3n3n1/10,
    A=−3n0/5+54n1n6/25.                           (6)

全十括弧の恒等式、行列式 δ≠0 を検証した。これも全源が一軌道である。

## 5. 全支持削減の証明書

`code/classification_degenerate_quadric_q3_anil3_support.py` は、
上記の全枝への元Jacobi代入、全三次式、補集合の因子、
三つの全パラメータ座標変更を厳密に検証する。
永続記録は `CLASSIFICATION_DEGENERATE_QUADRIC_Q3_ANIL3_SUPPORT_20261107.json`。
有限契約 `jobs/degenerate_quadric_q3_anil3_support.json` の下で1.71秒、PASS。
三つの canonical 代表と全源の括弧・座標変更・入力hashを保存した。

全支持の帰着により、代表一点の帰属を証明すれば全源を覆える。
有限点のJacobi成功だけから全源を推定したわけではない。
分岐・補集合・既約性の証明は統合担当が独立読解し受理した。
全三源のPGL恒等式も、§7末尾の独立Fraction証明書で全自由係数のまま再構成された。


## 6. upper 全源は T39 に入る

全五十次元の線形化kernelを、実軌道・Euler鉛直・scale の47方向と
真の三方向に分ける。相対トーラス重みは2,3,4で、全周囲の重み範囲は−6から8。
全150係数における真正103次元切片で、独立100式の陰関数は重み付き次数8以下に尽くされる。
有限な全陰関数を元の全350式へ代入すると、残余idealは厳密に零となる。
従ってこの点の射影局所次元は25で、正則かつ通る成分は一意である。
これは有限jetからの外挿ではない。

さらに canonical 代表を正則基点 `(a,b,c,r,s)=(1,0,1,0,1)` へ移し、
Darboux frame と全perspective変更を施して忠実S61座標へ正規化する。
全150係数のEuler商再構成と全134原方程式を確認し、
T39の**全被約graph**の45式が全て零になることを検証した。
T40のgraphには入らない。従って代表点はT39に属し、§1の全軌道恒等式から
upper の全二パラメータ源もT39に含まれる。

永続記録は
`CLASSIFICATION_Q3_ANIL3_UPPER_WEIGHTED_20261107.json` と
`CLASSIFICATION_Q3_ANIL3_UPPER_GRAPH_20261107.json`。
後者は全PGL行列、61座標、45残差を保存する。
帰属そのものには、全graph定理・忠実正規化・全源一軌道化で十分であり、
前段の局所正則性はその帰属と整合する別の確認である。

## 7. I・II 全源は正則 rank三二次因子への実曲線で MPS に入る

IとIIには真の二方向があり、相対重みは共に1,2。
全周囲の最大正重みは順に6,9なので、同変陰関数はそれぞれその有限重み以下の多項式。
全ての陰関数を求め、全350式へ代入すると残余idealは両者とも零となる。
こうして得た **実際の二変数Poisson族**を P(τ0,τ1) と書く。
全括弧は `CLASSIFICATION_Q3_ANIL3_CZERO_NONZERO_I_WEIGHTED_20261107.json` と
`CLASSIFICATION_Q3_ANIL3_CZERO_NONZERO_II_WEIGHTED_20261107.json` に保存した。

次の二曲線を取る。

    I:  Πraw(t)=P(t,5t²/2),
    II: Πraw(t)=P(t,t²/48).                              (7)

全tの射影Jacobiを元350係数で直接検証した。
さらに `Π=Πraw−E∧div(Πraw)/5` として標準unimodular liftを取り直すと、
全十strict Jacobiも恒等的に零である。
これは有限陰関数の存在だけに依存する主張ではなく、保存されたliteral曲線の全係数恒等式でもある。

実退化式は Q_t C_t で、二次因子はそれぞれ

    I:  Q_t=144q−144t ar+48t bc+8t²c²−t⁴r²,
    II: Q_t=186624q+31104t ar−10368t bc+288t²c²−t⁴r².  (8)

全四次minorは恒等零、(a,b,c)の三次minorは非零定数なので、全tで rank Q_t=3。
全五つのcanonical cofactor恒等式 `dQ_t·Π=Q_t X_t` と、
全五Pfaffianが実五次式のgradientであることを厳密に検証した。
cofactor行列の特性多項式は

    I:  λ(λ²−4t²)(λ²−16t²),
    II: λ(λ²−t²)(λ²−4t²).                              (9)

従ってt≠0でcofactorは正則半単純であり、非零固有値の二対の絶対比は2。
t=0のC_tは§4で証明した既約三次式の非零定数倍である。
可約な斉次三次式全体は射影係数空間で閉なので、0のZariski開近傍でC_tは既約のまま。
この開近傍のt≠0は `DEGENERATE_QUADRIC_CUBIC_REGULAR_SUPPORT.md` の
**全rank三・正則cofactor型**に属し、MPS成分の有限閉和に含まれる。
閉性からt=0も含まれ、§4の全源一軌道化からI・IIの全点が含まれる。

永続証明書は `CLASSIFICATION_Q3_ANIL3_CURVES_20261107.json`。
全曲線、unimodular lift、Q_t,C_t、cofactor、特性多項式を保存した。
`code/classification_degenerate_quadric_q3_anil3_curves.py` は全tの元係数恒等式を2.36秒で確認した。
`CLASSIFICATION_Q3_ANIL3_SUPPORT_REVIEW_20261107.json` は統合担当による独立Fraction再構成（0.54秒、PASS）。
全三源の普遍PGL恒等式、upperの全150商・134式・45graph式、
両曲線の全350 Jacobi・divergence gauge・strict Jacobi・全Q-normal/Pfaffian gradient・
全rank minor・特性多項式・t=0の既約性を再構成した。
帰属に必要な最小の算術入力はこのliteral二曲線、上の全源PGL恒等式、
upperの忠実graph認識であり、有限切片の全局所環解析を再検証する必要はない。
