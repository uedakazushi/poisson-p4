# dependent globalV全21族から六whole源への独立幾何監査

2026-11-07。reviewer `root`、producer `generic_exhaustion`。
**全必要係数、全例外、全逆変換の還元を独立確認しPASS。**
残六源のG上界又は全高重複度枝の完了は、本監査の結論ではない。
新しい係数標本、CAS、jetは使わない。

## 1. 二極端係数を別に導出した

normalssでa=t+1,b=t−1とする。
C項w²C∧dzに対し、full strict mixed Jacobiの寄与を項別に計算すると

`Jijz=(2b−a)w²(V∧C)ij, Jijw=0, Jizw=w³(bAC+kC)i`。

Euler wedgeからw³係数はx∧(bAC+kC)−(2b−a)Ax∧C。
z³のDv項も独立に同じ計算を行い、x∧(aADv−kDv)−(2a−b)Ax∧Dvを得た。
L,N,B,E等の全他項はnormal monomial z²w又はzw²側だけに入り、
weight二の全Π2及びhigher Jacobiはこのweight一方程式へ寄与しない。

normalnilでtrace tVの寄与は、同じ計算のa=b=t,k=0から
w³t*(x∧AC−Ax∧C)。zV∧dw及びαz²dzdwの寄与にはzがあり、この式を変えない。
必要な全C非零、normalssではDv非零を、実hyperplane因子という理由で保持した。

## 2. quotientによる全例外を含むlemma

x∧(pAC+kC)−qAx∧C=0、C非零、q非零を仮定する。
p非零ならW/Cの二次元quotientへ射影することでAC∈C-lineを得る。
p零ならim(kI−qA)⊂C-lineから同じrank-one修正表示を直接得る。
従ってA=cI+C⊗ell、trace零からell(C)=−3c、AC=−2cC。
原式に戻すとk=(2p+q)cとなる。
係数2p+qは反転していないので、その零補集合を失わない。
c非零はdouble semisimple、c零はA非零からrank一nilpotentである。

normalssの二qはt−3とt+3で同時に零にならず、全tで一方のlemmaが適用出来る。
初めにx=Cだけを代入してt=±2を除くような不十分な論法へ依存しない。
regularSS、mixed、Jordan3全族の不存在は全パラメータで成立する。

## 3. H/N全残余と正常化

H=diag(1,1,−2)では、C distinguished lineならk=3t−5、Dv distinguished lineならk=−3t−5。
Cのdouble eigenspaceはt=3,k=−2だけ、Dvのdouble eigenspaceはt=−3,k=−2だけ。
それぞれ他方の式がk=−14を要求するので成立せず、全場合でt0,k−5のみとなる。
N=E01ではc0からk0、tは全Cのまま、t=±3の拡大kernelも保持する。

normalnilのt0は既存別作者の全strict係数/trace-zero reviewを読み、Hα5/2とNα0を確認した。
t非零ではp=q=t,k0のlemmaからc0、従ってV=Nだけ。
normalnilのactual作用(T,V,α)→(σT,ηV,ηα)とbase nilpotent gradingを組み合わせ、
T=V=Nへ直した後もαの零/非零を0/1へ直せること、その全inverseを確認した。
α非零をt0の式の誤適用で捨てないというproducerの明示訂正は正しい。

これでnormalss H固定・N全t線、normalnil四固定の全六whole源を得る。
全symbol次元B=7,7(固定t6),6,5,7,6、plane六と射影化一の復元は、
元のdependent全Jordan正常形の全source数を戻したものであり、anchor軌道数だけを数えていない。
全30/9の残六源とparameter wallは次の独立の有限上界義務として保持する。
