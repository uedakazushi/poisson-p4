# rank零plane・全二十一低parameter源の全点被覆

2026-11-08、PRIMARY E、純全source証明。旧全350×30/9／全複素壁／全actual shear／全compatibilityの独立receiptを再利用する。新しいM1、profile、点、jet、GB又は全旧suiteを実行しない。独立純読解は後続receiptに委ねる。

**結論。** `PLANE_RANK_ZERO_GLOBAL_VECTOR_LOW_PARAMETER_WEIGHT_BOUND.md` に列挙された全二十一whole源について、rank四・被約既約actual quinticを持つ全点は、既essential三変数八23成分の閉和に含まれる。多くの源は既約quintic枝で空である。全rank零plane又はgeneric二／三parameter源を一括完了としない。

## 1. 全sourceと、元kernelによる不変z因子

planeはz=w=0、source30の順は十八base-pair項に続きD,B,E,C、bad Cはrows27..29である。全normalssのP0 z-row、及びnormalnil b零／非零のP0 z-rowはいずれもz倍。FULL9P2はbase-baseだけでnormal columnsを持たない。従って保存**全**K1上でCが零なら、全actualDにz因子があり指定既約枝は空である。SSではD零によるw因子も使えるが、この表には不要である。

既 `PLANE_GLOBAL_LOWPARAMETER_CLOSURE_SEPARATE_REVIEW.md` は全二十一源・全QQ parameterと全allowed wallでfull source/kernelを独立に再構成済みである。その全C rowを読むと、次の例外以外は全て零である。

* SS_COMMON_K0の固定fullK1。
* NIL_B0_TSQUARE_K0の固定fullK1。
* NIL_B0_COMMON_K0の固定fullK1。
* NIL_B0_MIX_TN及びNIL_B0_MIX_TN_VHの各tau=−8。
* normalnil b非零NDIAGのrho=−2/3。

最後のb非零例外は `PLANE_NORMAL_NIL_AFFINE_WHOLE_ACTUAL_FACTOR.md` の全rho／highm theoremにより既約枝で空。そのrho−1/3高m sourceも全C零である。generic rational coefficient fieldから壁内へmembershipを外挿せず、保存された全exceptional kernelそのものの零を使う。全minor／kernelのfactorは旧独立receiptが尽くしている。

## 2. T-square K零は全reduced compatibility上でz因子

旧 `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TSQUARE_20261107.json` は、全weight二compatibilityの被約台がu5=0、−u2+3u3+u4=0という単一linear spaceであることをexact generator identitiesから証明した。nilpotent厚みをunit idealとは呼ばない。

この代入後のfullP1 normal z-rowはu4*w*z∂a∧∂zだけで、全bad Cは零。他のz-rowもP0を含めz倍であり、全weight二kernel三係数を戻してもnormal columnは増えない。従って全actual Poisson点のDはz因子を持ち、このwhole源の既約quintic枝は空である。旧weight三非零式やheight上界を解き尽くす必要は無い。

## 3. SS/nil common K零の全源は本質三変数

旧fullweight二compatibilityの被約台は、SS commonで四、nil commonで三のlinear spacesの和である。全zero及び全intersectionを含む。旧独立proofは全actual finite shear inverse、全inhomogeneous P2 particularと全kernel三係数を戻してこの和へ移せることを確認した。

保存された全七branchのfullP1を読むと、一つも例外なく

```
P_cz=P_cw=P_zw=0,
P_ac=c L1(z,w)+Q1(z,w),
P_bc=c L2(z,w)+Q2(z,w),
P_az,P_aw,P_bz,P_bw はc,z,wだけに依存する
```

となる。ここではP0＋P1＋P2の全tensorをPと書いた。L1,L2は全binary linear、Q1,Q2は全P2由来binary quadraticで、zeroも全て許す。P_abはa,bを含む全式及び全三magnetic係数のまま残す。

実際、P0はSSでT=c∂a,V=c∂b,K0、nilでも同T,V,K0であり、上の形を持つ。全七branchのP1 ac/bcはc倍、c-normal rowsは零。P2は全base-base norm二なのでQ1,Q2を加えるだけである。

z一chartのintrinsic bivectorを(a,b,c,w)に取る。Euler補正を全て含めても

```
pi_cw=0,
pi_ac,pi_bc,pi_aw,pi_bw はc,wだけに依存する。
```

従ってactual signed quinticのchart表示は、向きの固定signを除き

```
D|z=1 = pi_ac*pi_bw−pi_aw*pi_bc
```

であり、a,bに依存しない。pi_ab及び全三magnetic項はPfaffianに寄与しない。dense chartでのpolynomial恒等式をhomogenizeして、全D=F(c,z,w)を得る。これはJacobiの未使用higher条件を零と仮定したものではなく、全必要sourceのtensor形だけからの結論である。

D零は指定枝で空。非零既約quinticなら本質二以下はC上のbinary factorizationに反するので本質三となる。後続の全curve/source/frame／六magneticを含むessential三変数八23成分のpointwise theoremへ送る。旧G像24上界だけからEを推論しない。

## 4. 二つのmixed tau−8の全sourceはc因子

T=N=b∂a、V=H+εN、H=diag(1,1,−2)、ε=1又は0、K=−8z²とする。旧全K1は八、J2はinjective rank九である。全六actual shearに対し、canonical normal blocksは

```
B=C=0, D=(kI−DV)s, E=(kI−DV)t.
```

kI−DV=−8I−(H+εN)はeigenvalues−9,−9,−6を持ち、非零定数det−486。従って全六D/E方向をactual finite inverseで消せる。残fullK1は、旧complete kernelの列3及び7の全二係数r,mとなる。

```
K3=−bw∂a∧∂b−(cw/3)∂a∧∂c+zw∂a∧∂z,
K7=(9aw+4εbw)∂a∧∂b+(4εcw/3)∂a∧∂c
       −3cw∂b∧∂c+9zw∂b∧∂z+w²∂a∧∂z.
```

全weight二はなおFULL9のまま保ち、finite shearのP2 translationを落とさない。Π0+rK3+mK7は**rawのまま**actual c→λcで固定される。c-rowは全てc倍、他coefficientsにcが無いからである。

実Poisson点のP2が存在すれば、旧J2のinjectivityで一意。c-scaleは全九P2空間と同じ方程式を保つので、一意P2自身も全λに固定される。P2_ac及びP2_bcはλでscaleするため零となり、残るP2_abの全binary quadraticを自由に保つ。よってfulltensorのc-rowは全てc倍であり、hyperplane c零がinvariant、actualDにはc因子がある。

これはP1の全r,mでcompatibleと仮定していない。実解が存在する点でのみ一意性とactual stabilizerを使い、その全sourceを可約又は零Dへ送る。global divergenceからstabilizerの延長を推論していない。r/m零、全P2及び有限normalizer逆を含む。

## 5. 全低parameter結論と、別の未処理lane

b非零は全rho/highmを含め空、generic／全壁でC零のsourceは空、T-square K零も空、二つのmixed tau−8も空。残るSS/nil common K零の全七branchは全essential三八成分へ入る。したがって宣言された二十一whole低parameter源の全点が既知閉和へ戻る。

normalss generic MIX三parameter、MIXH/MIXN/NILREG二parameter、normalnil b零Cartan三parameterとgeneric mixed二parameterはこの二十一packetの外であり、別source義務を保つ。scalar、co-Higgs k一、Π0零も本稿へ加えない。成分の新加算は無く、旧全report／kernel／manifestのbytesを変更しない。
