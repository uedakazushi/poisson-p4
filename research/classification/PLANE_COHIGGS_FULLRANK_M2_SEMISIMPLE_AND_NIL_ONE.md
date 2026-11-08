# Co-Higgs rank M三：全semisimple空性とnil tau一の全既知被覆

2026-11-08、PRIMARY E、log_manuscript。保存済み全kernelと実shearに基づくproducer純証明。
新operator、rank、minor、point、CAS又はjetを実行しない。
独立受理前のためglobal／1054主稿／frontier／Gitを変更しない。

## 1. 使用する全source

[元全rank三入口](PLANE_RANK_ZERO_FULLRANK_WEIGHT_BOUND.md)及び
[別作者の全係数review](PLANE_RANK_ZERO_FULLRANK_REVIEW.md)を使う。
U、M、T、Kの全actual normalizerと全逆を戻す。
M=[[b,c],[a,−b]]、n=(z,w)、En=z∂z+w∂wとし、全P0は

~~~
P0=T∧En+U∧((bz+cw)∂z+(az−bw)∂w)+K∂z∧∂w.
~~~

semisimple: U=∂b、T=2tau(a∂a−c∂c)、K=(1−2tau)zw。
nilpotent: U=∂a、T=tau(c∂b−2b∂a)、K=(tau−1/2)z²。
全weight一三十、FULL九P2、全inhomogeneous forcingと元350式を保つ。

## 2. Semisimple全tauは既約五次空

旧全QQ(tau) operatorと全rational根certificateでは、tau=0,1/2以外は
kerD1六=actual shear六、kerD2零。
選択minor固有の追加根tau=1,3/10,1/4,±1/6,±3/2,−1/2も、
全complete kernel／shear六／D2単射で個別に戻している。
[全finite shear補題](PLANE_COHIGGS_SHEAR_EXACT_SOURCE_REDUCTION.md)により
この全源はP0自身へ戻り、actual五次

~~~
−4ac tau(az²−2bzw−cw²)
~~~

のproper因子から空。
tau零は[全二normal定理](PLANE_COHIGGS_FULLRANK_HIGH_M_WHOLE_RANK_TWO.md)により全rank二以下。

残るtau=1/2でもkerD1六=actual shear六なのでP1を全finite逆で消せる。
D2完全核は九方向順(ab:z²,zw,w²; ac:z²,zw,w²; bc:z²,zw,w²)で

~~~
(0,0,−1; 0,2,0; 1,0,0).
~~~

従って全P2はrQ、
Qab=−w²、Qac=2zw、Qbc=z²に正確に等しい。
P1消去で誘導された二次P2補正も、この全一線へ入れた後の式である。
元のP0の全非零normal括弧は

~~~
Paz=az, Paw=aw, Pbz=bz+cw, Pbw=az−bw,
Pcz=−cz, Pcw=−cw, Pzw=0.
~~~

q=az²−2bzw−cw²と置く。
元signed Euler-Pfaffianを展開すると正確に

~~~
D(P0+rQ)=2q(rzw−ac).
~~~

手での確認として、a/c外側Euler項は各−acq、
z/w外側項は各rzwq、b外側項は零である。
qは非定数cubic、rzw−acは全rで非零の非定数quadric。
したがって全r、全zero境界で被約既約quinticにならない。
全actual正規化とfinite shear逆を戻し、
rank M三・M0 semisimpleの全multiplicity源が指定既約五次で空となる。
可約Dを持つrank四点まで不存在とは主張しない。

## 3. Nilpotent tau一：全P1を消してもFULL九の二核を保持する

旧tau一certificateではkerD1六、actual shear六、kerD2二である。
全finite shear逆でP1零へ移すと、全P2は

~~~
P2=(r z²−6s zw)∂a∧∂b+s z²∂a∧∂c,
        r,s arbitrary.
~~~

これは旧完全D2核二列
(1,0,0;0,0,0;0,0,0)、
(0,−6,0;1,0,0;0,0,0)そのものである。
全r,sを保ち、slice又は一点に置換しない。
normalized全tensorの非零上三角括弧は

~~~
Pab=r z²−6s zw, Pac=s z²,
Paz=−bz+cw, Paw=az−3bw,
Pbz=cz, Pbw=cw, Pzw=z²/2.
~~~

それ以外は零。全P2は共通∂aかつa非依存なのでJac(P2)=0、
旧D2核によってこの全二parameter族はactual projective Poissonである。
元signed Euler-Pfaffianは正確に

~~~
D=c²(az²−2bzw−cw²)+(r c/2)z⁴−s((b/2)z⁴+3c z³w).
~~~

以下の認識ではこの五次だけからtensor同一性を推定しない。

## 4. 全二parameterに共通する第二Poisson plane

S'=(c=z=0)では、全c-row／z-rowが消えるのでS'はPoissonである。
restrictionは−3bw∂a∧∂w、
S'上のEuler三vector係数は正確に3b²w。
一般restriction rankは二、零cubicはnonreduced 2L+Mである。

元actual DのS'に沿う最低normal次数は全r,sで三、
c,zをnormal variablesとしたその全leading symbolは

~~~
G3=−w c²(2b z+c w).
~~~

b,wが一般に非零なのでG3は恒等零でなく、normal位数は正確に三。
constant double root cとprimitive degree一simple root2bz+cwが一般に異なる。
r,sはnormal四／五次だけを変え、Poisson plane、rank、cubic及びm三を変えない。
basepointの特別点を取り除いてsourceのwhole性を失っていない。
これらは各tensorが満たすgeneric-along-S'の全入口条件である。

ここで既受理
[PLANE_RANKTWO_M3_WHOLE_POINTWISE_COVERAGE](PLANE_RANKTWO_M3_WHOLE_POINTWISE_COVERAGE.md)
をS'に適用する。
その定理は全R、全normal-root、全higher、全zero境界を尽くしており、
特に本2L+M/m三moving-simple枝をGa9既23／新HA23又は空へ入れる。
従って実Dが被約既約の全(r,s)は既知閉成分の有限和に属する。
全元sourceへfinite shearとactual normalizer逆を戻しても、
この既知有限和はPGL不変・閉なので全点包含が保たれる。

どちらの23成分名かをこの証明で一意に識別する必要は無い。
source全点を既にcatalogueへ採録された閉族の有限和へ覆うという、より弱いE結論で十分である。
新component数の増分は零。

## 5. 残余

rank M三・nilpotent側のtau=−3/2,1/2,3/2は本稿では未認識。
他のtauは旧全shear-exact補題、tau零は新全rank二定理へ戻る。
rank M二M(p)非零のm二全壁、rank M一以下、scalar、zero-symbolも別義務である。
今回の純因子／第二plane橋を全rank零plane又は全Eの完成へ外挿しない。
