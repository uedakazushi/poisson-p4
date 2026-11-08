# Normalnil b非零と全low21 sourceのpointwise E：独立監査

2026-11-08、PRIMARY E、別作者 `exceptional_coverage`。**全source純証明 PASS**。
対象は `PLANE_NORMAL_NIL_AFFINE_WHOLE_ACTUAL_FACTOR.md` と
`PLANE_GLOBAL_LOWPARAMETER_POINTWISE_COVERAGE.md`。
旧complete30/9、全kernel、全複素壁、actual normalizerと全multiplicity補足を読み、
新しい因子／essential-threeの論証を別に再導出した。
次元上界をmembershipへ読み替えていない。

## 1. 必要な固定byte読戻しだけを一回実行

job `jobs/plane_low21_pointwise_fixed_readback_review.json` gate PASS後、
標準Fraction/RFと旧独立checkerの基本疎多項式だけで固定arrayを読戻した。
0.05831秒、最大RSS26,100 KiB、全PASS。
producer/CASをimportせず、M1/J2/Jacobi/Pfaffian、rank、minor、kernel、point又はprofileは新生成していない。
旧証明済みkernelとminorの意味を再使用し、今回新たに使うbad-rowと全branchの形を直接確認した。

元source30はbase-pairのz/w十八項、その後D,B,E,C各三項。
bad Cはrows27..29、full9 P2はbase-baseだけでnormal columnsを持たない。
全二十一generic源と保存された**全元例外root**のcomplete K1を読み、C非零のlociは正確に次だけだった。

```
SS_COMMON_K0, NIL_B0_TSQUARE_K0, NIL_B0_COMMON_K0;
NIL_B0_MIX_TN / NIL_B0_MIX_TN_VH のtau=-8;
NIL_B1_DIAGONAL のrho=-2/3.
```

ほかの全generic/rootはC零。
normalnilのrho=-1/3高m rowは旧m二外flagを不変に保ち、全七kernelのC零も確認した。
rho=-2/3の追加列は正確にrow10=-4、row14=-1、row28=1だけ、J2 rank九/minor-96。
全七common branchではP0、全P1、**inhomogeneous P2＋free三係数**を足して、
c-normalとzw-row零、必要なmixed rowsのa,b非依存を直接確認した。
Tsquareのreduced support代入後の全normal z-rowもz倍だった。
実行全文、入力hash、全loci readbackを結果JSONへ保存した。

## 2. Normalnil b非零の全rho／全高m入口

旧独立source還元はnormal shearとtensor scalarの実逆によりK=zwとし、
行列commutator `[T,V]=V` から非零VはJordan3又はJordan21へ尽くす。
全centralizerの有限unipotent逆と共鳴rowの非零scale逆は全higherへ作用する。
従ってR、NDIAG全rho、NB、NCという全四型を保持する。

全multiplicity補足は位数二を使わずこのnormalizationを再適用し、
高mの唯一の残りrho=-1/3/resonant row零をNDIAG(-1/3)へ戻す。
旧certificateのexcluded_multiplicity_higher flagを変更せず、同じ全source/kernelへ適用する。
scalarなV零symbolやnormalnil b零を、この四型へ無断に追加しない。

P0のz-rowはz倍、fullP2はnormal columnを持たない。
C零の全sourceは従ってactualDにz因子を持ち、被約既約quinticでは空。
全rho generic、八つの元wall（高mを含む）、R/NB/NCを上記readbackで保持したので、
残るのはNDIAG(-2/3)だけ。
未調査のcomplex rootは旧全factor証明で無く、generic kernelをwallへ単純代入してはいない。

## 3. NDIAG(-2/3)の全六shearと残一係数

complete K1のC projectionはC_bの一列だけなのでker Cは六次元。
actual六shearは全てC零でrank六、従ってその像に等しい。
normal blocksは `B=-t,D=-(DV)s,E=s-(DV)t,C=0`。
したがって `t=-B,s=E-(DV)B` は非零sourceparameterを一切除さない全域の定数逆。
全finite shearはnormal weightを上げるだけで、P2への作用は全九変数のpolynomial translation、linear part identityである。

残一列は `-4bw∂a wedge ∂b-cw∂a wedge ∂c+w²∂b wedge ∂z`。
T=Ebase/3-b∂bなので、`E wedge(En/3+m*w∂a)` を引くと全P0+P1は

```
-b∂b wedge En+bz∂a wedge ∂w+zw∂z wedge ∂w
+m[-3bw∂a wedge ∂b+zw∂a wedge ∂z+w²∂a wedge ∂w+w²∂b wedge ∂z].
```

これは全mでc成分もc依存係数も無い。全九P2をなおA_ab+B_ac+C_bcとして保持する。
Euler regaugeはprojective Jacobiと同じJ2を保ち、旧rank九から存在するP2は固定m上で一意。
actual c-scaleはP0+P1をrawのまま固定し、full9 P2空間も保つ。
一意性によりP2もc-scale不変となり、ac/bc二binary quadraticは零。
abの全binary quadraticは制限せず残す。

よって全raw liftはc-row零で、actual signed Euler-Pfaffianは符号を除きc*Pf4。
零又はproper linear factorである。
全mで解が存在するとは仮定せず、存在する点でだけ一意性を使う。
実shear逆後の因子も非零linear formなので、全normalnil b非零・全multiplicityの指定既約枝は空。

## 4. Low21の残る例外を全て処理する

Tsquare K零の旧compatibility idealはu5平方を含み、u5零後は
`L=-u2+3u3+u4` に対してu3*LとL平方の恒等式を持つ。
したがって被約台はu5=L=0という単一linear spaceであり、schemeの厚みをunit idealとは呼ばない。
代入後P1のnormal z-rowはu4*w*zのみ、fullP2はnormal columns無し。
全actual点のDにz因子があり、残higher Jacobiを解かずに全空性が従う。

二mixed tau=-8では `T=b∂a,V=diag(1,1,-2)+epsilon*N`、epsilon=0又は1。
`kI-DV=-8I-DV` はdet-486で実六D/E shearに定数逆がある。
残二kernel列K3/K7はproducerの全式の通りで、c-rowはc倍、他係数はc非依存。
全九P2へのfinite translationを残し、同じinjective J2とactual c-scaleを使うと
P2_ac=P2_bc=0、全tensorのc-rowがc倍となる。
cはPoisson hyperplaneなのでactualDにc因子がある。
r,m零を含む全実解に対する論証であり、残P1を無条件compatibleと呼ばない。

残SS/nil common K零のweight二被約compatibilityは全四＋三branch。
全zeroとintersection、有限shear逆、全particular P2とkernel三係数は旧独立closure proofが保持する。
今回全七branchの全式から確認した形は

```
Pcz=Pcw=Pzw=0;
Pac,Pbc, Paz,Paw,Pbz,Pbw はc,z,wだけに依存する。
```

Pabは全a,b依存と全三magnetic項を保持する。
z一chartのintrinsic bracketはEuler補正込みで
pi_cw=0、pi_ac/pi_bc/pi_aw/pi_bwはc,wだけに依存する。
従って四次Pfaffianではpi_ab項が消え、actualDはa,bに依存しない。
dense chartの恒等式をhomogenizeして `D=F(c,z,w)` を得る。
unused higher Jacobiの消失又はmagnetic項零を仮定していない。

D零は指定枝外、非零既約quinticで本質二以下ならC上binary factorizationに反する。
ゆえに本質三で、既受理の全essential-three八23成分のpointwise theoremを適用できる。
全curve/source/frame／六magneticを含む同定が入力であり、旧G上界24から帰属を推論しない。

## 5. 正確な全点結論

宣言された全二十一lowparameter sourceのrank四・被約既約actual quintic全点は、
空な因子枝又は既essential-three八成分の閉和で被覆される。
normalnil b非零は全rho、高m、共鳴row零／非零、全higherを含め空。
旧all_passed=falseや旧G上界、旧m二外flagを成功へ書き換えない。

この二十一packet外のnormalss generic二／三parameter、normalnil b零Cartan／mixed一般源、
scalar、co-Higgs、Pi0零は別義務として保持する。
新component、全rankzero-plane E完了又は全Poisson分類を本稿から主張しない。
producer／旧証明書／global／main／Gitは変更していない。
