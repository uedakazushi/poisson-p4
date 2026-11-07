# F=x⁵+x²y³+y⁴zの全23次元成分：独立review

2026-11-07。E_reduced_nonnormal。**全source・周囲空間・真正支持・prior103相違 PASS**。
明示五次F=x⁵+x²y³+y⁴zから得る全cotangent sourceのPGL₅軌道閉包は、
projective dimension23の全Poisson既約成分であり、凍結した既存103成分全てと異なる。
成分数更新はintegratorに委ねる。他のfree quintic代表の尽くし又は全分類完成は主張しない。

## 1. 絶対既約性と全Saito source

Fはzに関して一次、leading coefficientはy⁴、constant termのy=0制限はx⁵。
したがってC[x,y]上primitive、C(x,y)[z]上既約であり、GaussによりC上絶対既約。
Fxのx⁴、Fyのy³z、Fz=y⁴を比べれば三偏微分は独立で、本質変数はちょうど三個である。
文献の例としての同定はproducerが記録しているが、この定義・既約性・free性の証明は引用に依存しない。

全quadratic monomial x²,xy,xz,y²,yz,z²を三成分へ置いた18方向と、
degree六の全28 monomialからgradient条件を別のFraction算術で構成した。
全28×18 entry、二kernel vectorの全行代入と独立性、独立に選んだnonzero16minorを確認した。
完全核は

`U=(0,−y²/4,3x²/4+yz)`,

`V=(−9y²/80,3x²/16−yz/4,9xy/40+z²)`。

独立determinant計算は`det(E_base,U,V)=−9F/64`。
既受理の代数的Saito criterionにより、これは次数二部分だけでなく**全AR(F)の自由基底**である。

全四frame係数Gと任意の六magnetic coefficient Rを残し、

`Π=U'∧∂u+V'∧∂v+[−div(V')u+div(U')v+R]∂u∧∂v`

の全十strict Jacobi、五divergence、五Casimir contractionを直接再展開して零と確認した。
全五raw Pfaffian及びactualDも照合し、`D=9 det(G)F/64`。
detG≠0だけがrank四fixed-F sourceの必要な局所化であり、全六Rを保持する。

既受理canonical全源定理とuniform磁気除去定理により、表示源はfixed-F全点を覆い、
実shearとvertex GL₂で一つのactual23軌道へ戻る。
R=0を単にgeneric magnetic sampleとみなしていない。

## 2. 同じ整数点から全150/350を再構成

standard frame・R=0を取り、分母の最小公倍数80を独立に求めた。
その整数点はproducerの六非零括弧と全150係数で一致する。
actualDは900F。
全150 elementary quadratic liftから全350 projective Jacobi微分Dを再構成し、保存された全entryと一致した。
producerのlegacy engineはimportしていない。

全Euler二十五、実GL₅二十五、scalar一の51列Mも独立構成し、全entryを照合した。
Mのexact rankは48、nonzero48minorとD M=0を確認した。
全fixed-F sourceは既受理の実一軌道なので、これは全sourceのraw orbit/gauge/scalar下界48である。

真正cutはactual48方向に非零な48座標を零とする線形cut。
保存48基底が実際のM列そのものに等しいこと、48座標上のdeterminant非零、
normal R₀がこのcutに入りD R₀=0であることを確認した。
101 independent implicit標準列と、actual48基底・R₀を合わせた全150直和も再検証した。
standard101列をLaplaceで除いた49×49determinantから、全150detは

`−3475607405985792000000000000000000000 ≠0`。

選択101 Jacobi行・101 implicit列のexact inverseについて左右積がidentityであることを確認した。
したがってrank D≥101。
一方、actual48とnormal一の独立49方向がker Dに入るのでrank D≤101。
よってrank D=101、raw tangent dimension49、projective tangent dimension24が正確に成立する。

## 3. actual C*と真正一変数IFT

座標重みσ=(0,0,0,1,1)の実PGL作用は、R=0点Π₀をt倍する。
全Saito sourceに対しても、pushforwardをtで割るとR→tRとして同じ全sourceを保つ。
従って正規化した実C*作用はΠ₀を固定する。
仮のweight assignmentでなく、actual stabilizerから得る作用である。

各raw inputのrelative weightはσi+σj−Σσk ek−1、outputはΣquad σi−Σσk ek−2。
全350×150非零entryでこのintertwiningを照合した。
input weight範囲は−3,…,1、唯一normal R₀のweightは−1。
actual48基底、coordinate cut、101implicit列、選択Jacobi行は全てこの作用でstableである。

選択101式のimplicit Jacobianが可逆なので、真正cutの完成局所環では
implicit座標yはnormal nの形式級数y(n)として一意に決まる。
解の一意性はC*同変性を与え、weight wの座標には−k=wとなるnᵏ項だけが入る。
y(0)=0、linear項はD R₀=0から零なので、事前に次数二・三までしかあり得ない。
ただし、以下の非零二次が既に最弱支持結論を与えるため、三次又は全graphは計算しない。

## 4. 全ambient余核の非零二次でreduced支持を閉じる

独立にQ=全350 Jac(R₀)を展開し、保存した九非零coefficientを照合した。
選択101式を消す一意の二次補正T₂の全150係数を照合し、
全350 residual D T₂+Q、選択101行の零、weight−2も確認した。
残差は四行で非零。

疎十一項functional λについて、**Dの全150列**でλD=0を再検査した。
そのexact pairingは

`λQ=−57897/220448 ≠0`。

これは選んだimplicit列だけの余核ではなく、全周囲空間の余核である。
任意の二次補正へ替えても同じ非零値となる。
真正一変数IFT残式の初項は非零n²なので、そのreduced local supportは原点だけである。
単に一つのfirst-order方向が持ち上がらないという一般多変数の推測ではなく、
full150直和で得た**唯一normalの真の切片**を使っている。

48linear cutのKrull不等式からraw local dimension≤48。
actual orbit/gauge/scalarの下界48と一致し、projective local dimensionは23。
PGL₅は連結なのでこの点を通る各全空間成分は23軌道を含む。
同次元と成分極大性により、その軌道閉包は全Poisson成分である。

tangent dimensionは24であって23ではなく、schemeのsmooth性は主張しない。
正確なnilpotence exponent又はfull completed local ringも求めていない。
必要なreduced成分性が二次で閉じた時点で停止した。

## 5. immutable prior103との全比較

比較台帳はcommit `228039f9f7e13ec850d0a49b0d21636fddbab814`、
`research/COMPONENTS.json` blob `e04cec24c8aed04e88053a57bb4a17cb8c2c1268`。
raw SHA256及びcommit→blobをGit read-onlyで確認し、projectionの全103行の全使用fieldを元blobと照合した。
可変の現行台帳にbindしていない。

異なるgeneric rank又は全成分dimensionのものは、それだけで異なる。
同じ23次元のMPS二十九成分については、受理済みraw/completed sourceの全periodを記号のまま残し、
全edgeを零にした**actual source**を再構成した。
全五次はκ x0x1x2x3x4、κは非零period式であり、gradient五次minorはκ⁵。
Laurent period chartとκ≠0の稠密開は非空なので、各component上の真正SNC essential5点である。
任意のlogcanonical polynomialを、そのcomponentの点として仮定したものではない。

新orbit closure全体ではactualDの五偏微分のspan rank≤3。
これはD係数のclosed determinantal条件で、Dが零になる境界も含めて閉じる。
従ってessential5のactual点を含む上の各componentとは等しくなれない。
一つの代表がgenericでないという議論ではなく、**新閉包全体**の閉不変量を使う。

既存の同次元既約quintic二成分についても、受理済みactual tensorからEuler/Pf五次を再構成した。

- R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23：gradient five-minor `−550376570880000000000`。
- R4_NIL22_LONG_IRREDUCIBLE_QUINTIC_23：gradient five-minor `−28476562500`。

両方とも非零であり、同じclosed essential条件で区別する。

残る同次元の非MPS成分については、immutable台帳の**既受理generic実D factorization**を入力に使った。
所定次数の積写像のprojective proper像はclosedであり、generic Dが可約/非被約なら全componentはその可約quintic閉錐に入る。
新しいactual900Fは絶対既約なので同一componentにはならない。
これら旧factorization証明を新たに再計算したとはしない。
全103個の理由は新review reportに個別保存した。
成分境界同士の交わりを排除する主張ではない。

## 6. 固定された独立証明書

finite gate後0.739秒、最大RSS23,608KiBで全exact reconstructionが終了した。
producer、legacy engine、CASをimportせず、既存の独立Fraction helpersを再利用した。
全入力は公開pathとSHAで保存し、旧artifactを変更していない。

|artifact|SHA256|
|---|---|
|code/classification_first_essential_quintic_review.py|4b6cf9b11124a4def824935c2769795ed955b03e88e02284bea310da2f278144|
|reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_REVIEW_20261107.json|c757b66b10196673d9b40370642d3e991e73309bd5f26d0e3a57e8b01c8f1a34|
|reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_REVIEW_JOB_20261107.json|9bfea029939281db3f8c3d2514ae97076824ed049cbd957df082bc2644a9bd05|

これは一つの明示Fの全source成分の証明である。
他のfree quinticの有限代表リスト・全代表の認識は次のE義務として残す。
成分数・global台帳・原稿・Gitはintegratorのみが更新する。
