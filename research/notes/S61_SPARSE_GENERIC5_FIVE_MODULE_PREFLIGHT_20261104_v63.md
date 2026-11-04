# 五生成加群による42個のunit座標targetと、global構成のpreflight

2026-11-04。実K5源上の次数四以下normal idealの次に、五生成加群からより小さい有限targetを取る。既存Fitting formulasと旧一個のdelta e0 targetは再発見として明記し、今回のactual coefficient solverや尽くしは未実行である。

## 1. 現行の閉残余と、旧結果との比較

R=K[n42]、A=E5×B、E5はR/m上の五次元有限étale代数、Bは原点を避ける閉埋込みR/I_Bとする。これらはNORMAL42 primaryと独立reviewの受理済み結論である。AはR上高々五生成であり、K上有限とは限らない。

旧NORMAL_SCALAR_FITTINGは既にFitt0(A)=m^5∩I_B、Fitti(A)=m^(5-i)(i=1,...,4)を与える。旧SCHUR_SQUARE_VECTORは、controlled localizationで一個のdelta e0所属と尽くしの同値も与える。本稿はそれらを新しい証明発見と呼ばず、actual42-normal baseのより低次数の代替targetと、source-bound構成の規模を整理する。

## 2. monic五次模型をnormal determinant無しで作る

任意のnormal依存行列のdeterminantを原点で非零だからと反転してはいけない。原点近傍のh-power basisやHensel逆行列だけではglobal Bを落とす恐れがある。次の構成では逆元をparameter係数体だけに限る。

既知planeのraw座標を(a,b,sigma)とし、係数体を忠実平坦にK(tau)へ延ばす。h=a+tau*b+tau²*sigmaとすれば、相異なる既知五点のh値の差は非零のtau多項式である。従ってhはE5のprimitive elementで、1,h,...,h^4がE5の基底になる。E5の一般分解として使うのはdegree3 toric étale代数＋K＋Kであり、さらに二次因子へ割れることは仮定しない。

plane idealは五二次式f5で生成され、その最高idealのHilbert関数は(1,3,1)。長さ5の一致からf5はfiltered H-basisでもあり、次数Dのideal元にはcofactor次数≤D−2を取れる。従って

```
y_l-p_l(h)=sum_i B_li(z)f_i(z), deg p_l<=4, deg B_li<=2,
g(h)=sum_i B_gi(z)f_i(z), gmonic5, deg B_gi<=3
```

をparameter-onlyな有限計算で構成できる。Macaulayを使う場合は35×55、56×105のfull-row-rank augmented行列だけで足り、その逆処理はK(tau)上である。normal変数を含むdeterminantを反転しない。

## 3. 優先するGram-freeな小さいplane構成

既存の五汎関数selectorをparameter行変換して、plane生成式を

```
G0=1+sigma-2sigma²+Q0(a,b),
G1=a+Q1(a,b), G2=b+Q2(a,b), G3=a*sigma, G4=b*sigma
```

へ正規化できる。三角形のtrace零と三つの異なる射影方向によりG3/G4のpure binary quadratic部分は零である。Q0,Q1,Q2の(a²,ab,b²)係数の3×3行列Mは非零detを持つ。canonical triangleではdet=3で、GL2変換により一般の場合も非零である。

toric factorの基底(1_T,a,b)ではbinary積ベクトルは−M⁻¹(1_T,a,b)。これからa,bの二つの3×3乗法行列を得て、h_T=a+tau*bの3次特性多項式g3を取る。3×3 Krylov基底(1_T,h_T,h_T²)は可逆。elliptic二点のh値は−tau²/2、tau²なので

```
g(h)=g3(h)*(h+tau²/2)*(h-tau²)
```

である。p_b,p_sigmaはtoric degree2表示とell二点のChinese interpolationで得られる。必要な逆元は3×3 determinant二つ、g3(−tau²/2)、g3(tau²)、tau等のparameterだけで、すべて非零はprimitive hから従う。

座標復元とhの定義も多項式として一致させる。bとsigmaのcofactorだけを先に選び、
`p_a=h-tau*p_b-tau²*p_sigma`、`B_ai=-tau*B_bi-tau²*B_sigma_i`
と定める。これで三つの補間恒等式は同時に成り立ち、cofactorにも`B_a+tau*B_b+tau²*B_sigma=0`が厳密に成立する。後のnormal correctionにも`S_a+tau*S_b+tau²*S_sigma=0`が成立するため、hの定義を追加のrelationにする必要はない。三つのcofactorを独立に選んだだけでは、この最後の一致を仮定してはいけない。

cofactorはgraded lex(a>b>sigma)の六border rewrite a²,ab,b²,a*sigma,b*sigma,sigma³で追跡できる。最初の三つはM⁻¹、次二つはG3/G4そのもの。最後は

```
sigma³-(sigma+sigma²)/2
=-(sigma/2)G0+((q00*a+q01*b)/2)G3+(q02*b/2)G4
```

(Q0=q00a²+q01ab+q02b²)である。標準単項式は1,a,b,sigma,sigma²の五個。五点代数の長さと一致するので六rewriteは全idealを生成し、degree≤4/5 targetの除算は高々35/56単項式の範囲で終了する。新しいGB solverを呼ぶ必要はない。

## 4. global normal correctionと、自由五次代数

cross rewrite z_h n=R_h(n)nを使う。R_hはaffine-linear normal行列として、既存のcanonical Lと実130法線因子からsource-boundに定義できる。選択五式はcrossを減らすとf_i(z)+K_i(n)、K_i(0)=0、degree≤2になる。

K_i=k_i(n)nと因数化すれば、z-wordを掛けたK_iのnormal reductionはk_i R_reverseword nとなる。これは全A内の正確な恒等式であり、normal行列のscalar係数を展開せず長さ≤3の行列積として保存できる。cofactor identitiesから

```
S_l=-Red(sum_i B_li K_i), degS_l<=4,
H=-Red(sum_i B_gi K_i), degH<=5,
y_l=p_l(h)+S_l(n), G(h,n)=g(h)-H(n)=0
```

を得る。これはnormal0 fiberやformal completionだけの関係ではなく、全Aでの多項式恒等式である。

従ってC=R[h]/GはR自由rank5でAへ全射する。元の選択五式＋130残差を代入し、Gで割った135個のremainder r_iを取るとA=C/(r_i)。逆写像では元の座標をp_l(h)+S_l(n)へ送り、前節のcofactor選択によりa+tau*b+tau²*sigma=hが恒等的に成立する。このため余分なh-definition relationはなく、双方の写像は互いに逆である。全675列rem_G(h^k r_i), k=0,...,4がactualな5×675 R-module presentation Dを与える。monic Gのcompanionの多項式としてだけ乗法行列を作るため、任意のlift行列の可換性を仮定していない。im DはC-idealであり、h倍で閉じる。

normal次数を保ったまま見積もると、130残差の650列はdegree≤6、選択五式の25列はdegree≤10。理由はh^mのG剰余に現れるHの次数が高々floor(m/5)で、残差はh-degree≤4かつnormal-degree≤5、選択五式のh-degree≤8/normal-degree≤8だからである。675行列の実係数はまだ構成していない。

## 5. 42個のunit座標targetが完全な尽くしcriterion

e0=(1,0,0,0,0)はCの単位1へ写る。次は正確に同値である。

```
B=0
iff for all j=0,...,41, there exists U_j in R^675 with D U_j=n_j e0.
```

右辺ならAでn_j*1=0、従ってmA=0となりA=A/mA=E5。逆はA=E5なら全n_jが零なので、完全presentation Dから各cofactorが存在する。210個のn_j e_k targetは不要である。im DはC-idealだから、unit targetをh^k倍すれば残りをすべて得る。

次の検証の完了条件は、42本の**有限多項式**U_jを得て全normal係数についてDU_j−n_j e0=0をQQで示すこと。一次jetや形式級数の可解性では足りない。cofactor次数の上界は現時点では証明していない。

### 同じformal tauで一個の線形targetへ符号化する

同じ独立変数tauを用いてell(tau,n)=Σ_(j=0)^41 tau^j*n_jと置くと、上記は一個の有限module target

```
exists U in K(tau)[n]^675 with D U=ell(tau,n)*e0
```

とも同値になる。元のAはtauと独立にK上で定義されている。K-vector spaceの基底を用いればA[tau]→A⊗_K K(tau)は単射である。従ってell=0ならtauの全係数n_jがAで零となり、A=A/mA=E5である。逆も明らかで、完全module presentationからcofactorが存在する。

この証明にnormal linear formの逆元や、B上のregularityの追加仮定は不要である。normal次数は1、tau次数は41。ただしこれは42個のtargetを一個のformal多項式へ符号化したもので、cofactorのtau次数・normal次数の上界を与えない。旧delta e0よりsolver costが必ず小さいとは主張しない。Dの表示をtauに依存させても、表すA自体は忠実な係数拡大なので議論は変わらない。

## 6. Fitting側の有限colon5と限界

E5はR加群として(R/m)^5、B=R/I_Bでm+I_B=Rなので

```
I5(D)=Fitt0(A)=m^5 I_B=m^5 intersect I_B,
Fitt_i(A)=m^(5-i), i=1,...,5.
```

旧formulaの直接のcorollaryとして

```
I_B=I5(D):m^5,
B=0 iff m^5 subset I5(D)
```

を得る。colonは不特定の飽和指数から有限指数5へ固定されるが、maximal-minor idealの実所属自体は未計算である。m^5の単項式は1370754個あり、42個のvector targetの方を次の実装候補にする。低いFitting idealsはBの有無によらず同じで、entry ideal=mだけでは尽くしは従わない。

## 7. associativity・五生成だけでは足りない反例

canonical plane f0,...,f4に対して、K[a,b,sigma,n]で

```
(f0-n,f1,f2,f3,f4, n*a,n*b,n*sigma,n²-n)
```

を取る。nはidempotentで、n0 fiberはE5、n1 fiberはa=b=sigma=0の一点。従ってAtoy=E5×K、非空のBtoyがn1に残る。crossはfull rank3、canonical plane HF(1,3,1)、R=K[n]上五生成、commutative/associativeをすべて満たす。

h=a+4sigmaを使うfree五次coverのh根は1,0,-1,-2,4。R-moduleは四個のR/(n)と一個のR/(n(n-1))で、diagonal presentationのFitt0はn^5(n-1)、高いFittingはnの通常冪になる。形式的局所解とglobal残余の違いを具体的に示すtoyであり、実βの数値標本ではない。

## 8. 次の一回の候補と実行しない範囲

次に提案する有限計算は、実K5の50個のplane係数と既存selectorから上のM3、二つのtoric乗法行列、formal tau Krylov、g5/p、六border cofactorだけを**field SLPのまま**構成・点検する一回である。これなら35/56 Gramさえ不要で、最大56個のplane単項式、3×3/2×2 parameter逆元だけを扱う。既存126 Gramのscalar expansionには触れない。

その後のD normal係数や所属探索には別予算が要る。一次normaljetだけでも210×675=141750個のparameter scalar slotがあり、これはparameter多項式項数の上界ではない。denseな42-normal degree10展開の一般上界は巨大であり、許可なく展開・Gaussian・GBを始めない。今回のscriptはsource binding、degree/support preflight、toy CRTの有限確認のみで、M3の実値、g/p、D、U_jはまだ計算しない。

## 9. 独立に小さく検証できる次のinterpolation計算

actual K5のM3係数を直ちにQQ五変数分数へ展開する上界は150000項と32MiBを保証しない。そこで次の許可候補は、canonical triangle自体を使う**三つの形式係数**alpha,beta,gammaだけの計算である。h=alpha*A+beta*B+gamma*sigmaとすると既知五値は

```
theta=(alpha,beta,-alpha-beta,-gamma/2,gamma).
```

g5(h)=Π(h-theta_i)、p_A,p_B,p_sigmaは明示Lagrange和で得られる。共通分母にはVandermonde D=Π_(i<j)(theta_i-theta_j)を使い、Dで払った三つのinterpolation identityとg5(h)を、定数係数のcanonical六border rewriteでcofactor付きに還元する。ここに新しいβ点・素数・Sp4不変式の次数探索はない。

三つのcleared targetはparameter次数10・plane次数≤4、g5 targetはparameter次数5・plane次数5。従って、各targetの係数support上界は66×35と21×56、各cofactorは66×10または21×20である。全target・cofactor・小さいg/p出力を含めた保守的上界を**30000項・有理係数512bit以下・RAW12MiB以下**とする。定数border rewriteは各targetで高々56個のplane単項式を処理し、6の冪までの係数分母増大を含めてもこのbit上界内に収まる。

height boundも具体的に取れる。2thetaの差の整数linear formのL1 normは高々6で、Dを払ったtargetのnumeric denominatorを2^11で揃えると整数係数は32bit以下。canonical M3の2×2 cofactorは8以下、逆行列の分母は3で、各border replacementのL1増大は16以下となる。56段後の共通numeric denominatorは2^11*6^56、cofactorを蓄積する係数のbit数は32+56 log2(16)+56 log2(6)+10<512で抑えられる。これは係数上界の証明であり実係数を計算した結果ではない。

actual sourceとの対応は、toric triangleを(a,b)=(xA+yB,zA+wB)とするとalpha=x+tau*z、beta=y+tau*w、gamma=tau²である。三角形のGalois permutationでLagrange式は一意なinterpolationとして降下し、必要なら六個のS3変換によるcofactorの平均を取れる。actual tableの3×3/Krylov/Chinese有限recipeとの対応も、ここで独立に確認できる。既存126 Gramは一切評価しない。

これは次の一回の**候補とsupport preflight**であり、今回その三形式係数計算を実行したわけではない。normal correctionや675行列のscalar係数に進む予算もまだ与えていない。

## 確定した保存範囲

今回のsource/theorem preflightとtoy確認は0.4019秒で完了した。source-bound recipe、有限colon5、42個と一個のformal linear target、既存成果との違い、反例模型を保存した。追加の一個targetとinterpolation予算注記は数学的・metadata追記だけで、actual field solveやnormal係数計算を再実行していない。現在の成果は新しい所属証明そのものではなく、その正確な完了criterionと、scalar Gram展開を避ける次の検証単位である。
