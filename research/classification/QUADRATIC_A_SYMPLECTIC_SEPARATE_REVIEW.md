# quadratic_A全七源のsymplectic MPS被覆・独立review

2026-11-07。主義務E_multiple_hyperplane。producerと別作者による全source幾何及び独立exact算術。**PASS**。
**quadratic_Aのλと全六磁気係数を持つ全源、そのPGL飽和及び閉包は、既知MPS40の有限閉和に含まれる。**
新成分を加えない。停止中の十八変数局所ideal、全形式変形又は全symplectic源の分類を結論しない。

## 1. 全七源と実可逆normalizer

α=λ−1/2、β=−λ−1/2、I=m0+m4、
R=m0+m1x+m2y+m3x²+m4xy+m5y²とする。
元の全源は

`{q1,q2}=R+βyq1−αxq2`,

`{q1,x}=αx², {q1,y}=1+βxy, {q2,x}=1+αxy, {q2,y}=βy², {x,y}=0`。

独立checkerは公開全source reportのcase7と、λ及び全六mag parameterを記号のまま全係数で照合した。
元の全四strict Jacobiも展開した。

αβI≠0で

`q1'=(q1+m2/β+m5y/(2β))/I`,

`q2'=(q2−m1/α−m3x/(2α)+m4y)/I`, `x'=x, y'=y`

及びtensor scalar Iを用いると、全六括弧が同じλでR=1へ戻る。
逆は

`q1=I q1'−m2/β−m5y/(2β)`,

`q2=I q2'+m1/α+m3x/(2α)−m4y`, `x=x', y=y'`。

両inverse composition、Jacobian determinant I⁻²、全六pushforward恒等式を独立QQ(λ) Laurent算術で確認した。
α+β=−1を正確に保つため、Iを独立変数にしてm0=I−m4と書き、六magを失っていない。
この変換はaffine linearであり、hを固定する実PGL₅変換へ斉次化出来る。

次に(U1,U2,V1,V2)=(q1',x,y,y−q2')で定数部分は
J=∂U1∧∂V1+∂U2∧∂V2となる。
全sourceはR=1について正確に

`{U1,U2}=αU2²`, `{U1,V1}=1+βU2V1`,

`{U1,V2}=−U2V1−βU1V1−αU2V2`, `{U2,V1}=0`,

`{U2,V2}=1+αU2V1`, `{V1,V2}=βV1²`                       (1)

へ移る。この全standardizationも六括弧で照合した。

## 2. 真正全parameter親族とSNC定理の適用

最弱証明に必要な親は、自由な三parameter A,B,cを持つ次の実Poisson源だけである。

`{u1,u2}=−c u1u2`, `{u1,v1}=1+A u1v1+B u2v2`,

`{u1,v2}=c u1v2`, `{u2,v1}=(B−c)u2v1`,

`{u2,v2}=1+B u2v2`, `{v1,v2}=(B−c)v1v2`。                 (2)

全四JacobiはA,B,cを自由に残して恒等零。
Pfaffianは`−(1+B u2v2)(1+A u1v1+B u2v2)`。
AB≠0では斉次二因子は相異なるrank三・rank五quadricであり、C上既約である。
ただし全H+2Q分類を依存として仮定する必要はない。

独立checkerは重み(h,u1,u2,v1,v2)=(3,0,1,0,1)を全項へ適用した。
全項重みは0,2,4,6で負のものは無い。
正重みの定数項及び交換項だけが消え、極限はlogcanonical、実五次式は

`−AB*h*u1*u2*v1*v2`。

AB≠0だから五独立超平面の真正SNC log-symplectic点である。
cの値、modular eigenvalue ratio又はgeneric periodを追加仮定していない。

これは既受理 `TWO_QUADRICS_SUPPORT.md` §3Bの全源と同一である。
座標(x,y,z,w)=(u1,v1,u2,v2)で

`C=D=1,E=0,F=B,p=−c,q=c,r=c−B,s=B−c,v1=0,v2=−B`

という**全六括弧の同一性**も独立に確認した。
connected PGL₅は全空間の各irreducible componentを保存するため、元の親の点を通る各成分はそのSNC極限も含む。
既受理MPS Theorem1.4の全period・全SNC成分のscopeにより、親は既知MPS40の閉和に入る。
SNC点での局所成分一意性や特定のMPS番号は使わない。

## 3. literal曲線、全λ端点、全分母補集合

(2)にA=−ε(2ε+1)、B=−ε(ε+1)、c=εαをそのまま代入し、
ε≠0で

`U1=u1+u2/ε,U2=u2,V1=v1,V2=v2−v1/ε`                  (3)

を作用させる。U側G、V側G⁻ᵗの実変換であり、determinant一、Jを正確に保存する。
逆はu1=U1−U2/ε、v2=V2+V1/ε、u2=U2、v1=V1。
全pushforwardを元の真正親から再展開すると

`{U1,U2}=αU2²−εαU1U2`,

`{U1,V1}=1−ε(2ε+1)U1V1−ε(ε+1)U2V2+βU2V1`,

`{U1,V2}=εαU1V2+(2ε−β)U1V1−αU2V2−U2V1`,

`{U2,V1}=ε(β−ε)U2V1`, `{U2,V2}=1−ε(ε+1)U2V2+αU2V1`,

`{V1,V2}=(β−ε)V1(V1+εV2)`。                          (4)

全負冪は相殺し、ε,λの文字通りの多項式族となる。
ε=0は**全λ**について(1)と一致する。
独立checkerは全四strict Jacobi、hによる全150斉次係数及び全350 projective Jacobiを検査した。
実Pfaffianも元の親、変換後の族、公開式の三者で一致する。

AB≠0の条件はε∉{0,−1,−1/2}だけであり、λのwallは無い。
この稠密curve chartの全点は§2の実SNC定理で既知閉和に入る。
有限MPS閉和の閉性によりε=0も入る。
他の二つのε例外も、curveがregularなので閉性で含まれるが、端点の認識に別計算は不要である。

発見時のpole説明も独立に確認した。
一般のspectral係数r,sと真正period cを保って(3)を作用させると、r,s係数の最大poleは二、c係数は一。
r=ερ0+ε²ρ1、s=r+ε、c=c0+εc1とすると、正しい全四unknown系のrankは四で、
(ρ0,ρ1,c0,c1)=(1,2,0,λ−1/2)が全十八行を正確に満たす。
二次spectral項は端点に寄与し、それより高いspectral項又は二次以上period項は端点に寄与しない。
ただし最終帰属はliteral族(4)だけで証明出来、全contractionの完全性を要しない。

## 4. 稠密七源から全補集合とPGL飽和へ

元のλ,m0,…,m5のaffine parameter space A⁷は既約である。
全sourceの射はpolynomialであり、定数括弧が非零なのでprojective tensor空間への射も全域で定義される。
αβI≠0は非空稠密principal open。
その像は§1の実座標同値と§3の全λ認識からMPS閉和へ入る。
既知閉和の逆像はA⁷のclosed subsetなので、**α=0、β=0、I=0を含む全七源**を含む。
MPS閉和はPGL₅不変かつclosedなので、全飽和及びその閉包も含む。
chosen magnetic値から全磁気を推論したものではなく、全六係数normalizerと稠密閉性を用いた全源の矢印である。

## 5. 補助五係数源と最弱依存の区別

公開source reportにある表示五係数源についても、独立に全Jacobi係数を展開した。
b=q+r、e=p+qに対するideal生成元は

`be, b(b−f), e(e−a)`。

保存されたconstant coefficient matrixの両方向span恒等式を再検証した。
従ってこの**表示源**のreduced supportはb=e=0、e=0/f=b、b=0/a=eの三parentへ尽くされる。
これは有用なexact代数証明だが、上の最弱帰属証明には不要である。
全generic semisimple源のweight完全性や全非semisimple親分類を新たに仮定せず、literal親(2)だけで閉じる。

旧三unknown負結果及び停止中の局所idealのreceiptは変更していない。
負結果を全source非所属又は新成分の根拠としていない。

## 6. 固定証明書

finite job gate PASS後、0.383秒、最大RSS13,824KiBで全算術が終了した。
producer及びCASをimportせず、既存の別作者Fraction(λ) field実装だけを再利用した。
全依存path/hashはjob及びreview reportに保存した。

|artifact|SHA256|
|---|---|
|code/classification_quadratic_a_symplectic_review.py|8d76f227966434e76db7ad63bfb80575846a7591d71c2bec37f920b118028b1a|
|reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_REVIEW_20261107.json|566de4208dda9962280f6069ca5b9b898f846f92687c61900bdd13343110535c|
|reports/CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_REVIEW_JOB_20261107.json|e9f5727178185b94192dab82762fe3e46d367c7904facfba295afb607017732a|
|producer CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SOURCE_20261107.json|d44e818148f74e4bc32aa176e967b0569fc17eb86f85afeeef809da6581af4c0|
|producer CLASSIFICATION_QUADRATIC_A_SYMPLECTIC_SPECTRAL_ARC_20261107.json|70a273c5849e2a3df8ca0fb4b5075702cac8e470c604ec2ac43d1d12ad7ad45b|
|producer QUADRATIC_A_SYMPLECTIC_MPS_LIMIT.md|2003b98e72a0355e4f1bf4100afe6b09060e7df87774cca40ab53c2aea793c9a|
