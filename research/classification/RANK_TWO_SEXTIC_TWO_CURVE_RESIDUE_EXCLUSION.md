# 次数六rank二支持の全空性：二種類の曲線上の留数定理

2026-11-08、PRIMARY E、rootによる新しい純証明。別作者読解前。
新CAS、射影点、operator、kernel、minor又はjetは無い。

**定理。** rank四Poisson tensor on P4のactual退化因子Dが被約既約quinticなら、
非正規支持にdegree六integral surface Sがあって、その一般restriction rankが二であることは無い。
全有限射影・node/cusp・split/non-split・small resolution・codimension二導手境界を含む。
旧G次元上界だけをE membershipへ読み替えた主張ではない。

## 1. 既受理の全source入口

`RANK_TWO_SEXTIC_SURFACE_SEPARATE_REVIEW.md` §§2–4を使う。
length/Bertini/Höringによるrank非依存の全入口から、Dの正規化は完全系Lのdegree五minimal scroll
S(1,1,3), S(1,2,2), S(0,1,4), S(0,2,3), S(0,0,5)の全五型である。
他の二次元非正規支持は無く、δ_S=1である。

005はline-vertex coneの像で、二次元非正規支持はplaneの有限和だけとなり、degree六のSは不可。
他四型ではsmooth threefold
Y=P_P1(O(e0)⊕O(e1)⊕O(e2))を使う。113/122はY自身が正規化、014/023は
zero-summand sectionだけを縮めるsmall resolution Y→Xである。
縮むdivisorは無い。Pic(Y)=ZL⊕ZF、K_Y=-3L+3Fで、divisorial conductorのstrict transformは

    C ~ 3L-3F。

cuspはC=2ΓとなりL係数三の奇数性に反する。
従ってgeneric nodeで、Cは一prime Γ又は二prime Γ1+Γ2。
一primeにはnon-split nodeのconnected global conductorも含む。

同じ旧全源定理により、正規化上のclosed inverse-residue一形式αは、
全codimension一でregular Poisson tensorとdualizing generatorのcontractionである。
conductorに沿って高々simple logarithmic poleを持ち、各primeのresidue λは非零定数である。
非零性はrestriction rank二とSへのgeneric separabilityから従う。
係数体をPoisson中心とは仮定しない。
small resolution上でも、exceptional curve外のαを有理一形式として取り、
全codimension一で同じpole divisorだけを持つ。新しいexceptional divisorは生じない。

## 2. 曲線への制限で使う補題

smooth projective threefold Y上のclosed rational一形式αが、prime divisor Γiに沿ってのみpoleを持ち、
genericにsimple/logarithmic、residue定数λiであるとする。
smooth complete curve Bが各Γiのgeneric smooth locusを横断し、
Γi相互の交差、singular locus及びその他のbad codimension二集合を避けるなら、
α|Bは他にpoleを持たず、その全residue和は

    Σ λi (B·Γi)=0。

理由は局所横断parameter zに対しα=λi dz/z+regularであり、
各交点のpullback residueがλiとなることと、compact curveのresidue theoremである。
接触multiplicityを許す場合も同じ式だが、以下は横断曲線だけを選べば足りる。

ここで必要な二種類の曲線はいずれも存在する。
第一は一般のL-hyperplane二つのcomplete intersection B_H。
small mapのvertexと全bad codimension二集合を像で避けて選べ、smooth curveとなる。
第二は一般のruling fiber P2⊂Y内の一般のline B_f。
bad集合の各curveがbase P1をdominateしても、一般fiberとの交わりは有限集合だけである。
finiteな点を避け、Γiのplane curveと横断するlineを選べる。
exceptional zero-sectionとも一つの点でしか交わらず、それも避けられる。
全curveがfiberであるbad集合は有限個のbase値だけを占め、一般fiberで避ける。
従ってsmall resolutionのcodimension二を無断で無視したのではない。

## 3. 一primeの全排除

C=ΓならL²Γ=12である。
B_Hへの制限と補題は12λ=0を与える。
λ≠0に反し、一prime/non-splitの全場合は空。
これは旧独立受理済みの排除をそのまま確認したものである。

## 4. 二primeに対する二つの独立した留数関係

C=Γ1+Γ2なら、nodeの二global branchはそれぞれSへbirationalに写る。
従って

    L²Γ1=L²Γ2=deg S=6。

B_Hに補題を用いると

    6λ1+6λ2=0、従って λ2=-λ1。                  (H)

一方Γi=a_i L+b_i Fと書くと、a_iはruling P2に制限したdivisorのdegreeであり整数である。
全conductor classから

    a1+a2=3。

B_fに同じ補題を用いると

    a1 λ1+a2 λ2=0。                             (f)

(H)を(f)へ入れると(a1-a2)λ1=0。
a1+a2が奇数三なので整数a1,a2は等しくない。
従ってλ1=λ2=0となり、rank二に必要な非零residueに反する。

特定の射影、mark、conductor sectionの因数分解又は一般位置のsourceを選んでいない。
|2L-4F|にprimeがあるかどうかさえ必要ではない。
旧残余023の全split-conductor proper projection locusにも同じ二関係が適用する。
一つのL²曲線だけの留数条件では不十分だった点を、ruling lineの独立条件で閉じた。

## 5. 全sourceへの復帰と射程

005、cusp、一prime、二primeを合わせて全degree六rank二の既約quintic枝は空。
smooth/non-normal支持S、全finite五section、small resolution、全codimension二導手境界は
上の旧全source入口と曲線の選び方に含まれる。
旧G非優越証明・witness・raw certificatesを変更しない。
新結論はそのproper residue必要条件を全点で強化するもので、旧算術の誤りを主張しない。
rank零degree六は別既受理排除に従う。全E、他degree、他factor又は総成分数は未主張。
新成分増分は零。別作者がこの全source適用と二curveの存在を読解するまでglobal採録しない。
