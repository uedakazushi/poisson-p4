# cubic scroll A25 全境界証明：別作者読解と一箇所の正確なscope補正

2026-11-07、主義務 G、reviewer `/root/log_manuscript`。
**新しいboundary証明は PASS。** 対象 producer は
`CUBIC_SCROLL_A25_WHOLE_NONDOMINANCE.md`、SHA
`f08455bf52022f7cb3932f9428801c204d969189bd3903ed66a6a6f2c6c999a5`。
producer と全入力 bytes は変更しない。新CAS、係数sample、minor、jetは実行していない。

## 1. 独立性と受理する正確な範囲

reviewer は、producerが入力とするscalar-connection、quadric-net geometry、net-Poisson reduction、
exceptional-restriction uniquenessの四つの先行noteの著者である。
従って本receiptはそれらを「完全別作者によって再認証した」とは呼ばない。
それら四noteは今回のproducer `plane_collision_review` が別作者として全読解し、
本reviewerはその先の新しい §§4–10 の全boundary/classification/count proofを別作者として再導出した。
root integratorによる全chain読解と組み合わせる、明示した相互読解の分担である。

受理する結論は、smooth cubic scrollの残る二normal-root型
`A=2C0+5f` のrank一、及び同Aのrank二pairについて、固定scrollのprojective Poisson像 dimension≤6、
全PGL5飽和像 dimension≤24である。
canonical-root排除と結合してsmooth cubic scroll rank二の全G非優越へ適用できる。
全点のE帰属、cubic cone、全G/全分類、新成分数は結論しない。

## 2. 全sourceとextensionの境界

先行scalar reductionはRの全multiplicitiesを保持する。
rank二のeta pole budgetは、detM order mとbeta pole order nをprimitive unitで比較して n−m。
R=C0+divdetMからetaはC0上simpleだけとなり、GysinとH0Omega1=0で消える。
rank一の第二行はnormal rootとして選んだものではなく、residueのprimitiveから得るglobal sectionである。
**そのdivisorial zeroを許す**ことをproducerが明記している点を確認した。

rank一のzero-residue primitiveについて、birational座標の仮定を避ける別証明も確認した。
proper rational EはH1(E,C)=0、divisor complementのH1はmeridiansで生成される。
全residue零により全period零、global analytic primitiveがあり、generic divisorで有限poleを積分すればmeromorphic。
codimension二Hartogsでproper E全体へ延長し、meromorphic functionはrationalとなる。
alphaのfiber成分sigma/Gとdfのfiber numerator比較はdetM=const*sigma/gammaを与える。

全Y=P(F)はsmooth P2-bundle。D'はH+2Bのrelative hyperplaneであって、その全-zero basepointではfull P2 fiberを持つ。
L上のF=O²+O(−2)とV=C0から全binary quartic kappaが残り、kappa≠0はDの既約性による。
四つのzeroの重複、無限遠、L外の追加full fibersを一切除外していない。
normalization-flat性も仮定していない。

ambient extensionを数える際はH0(Y,Lambda²TY(−E))=0が必要であり、producerはこの全域vanishingを使用している。
単なるformal first jetからglobal Piを推測していない。
pi_Eはkerdfとsigmaから一意、よってrow-pair parameterからprojective Piは高々一つである。

## 3. dual determinantと全relative hyperplane

全六係数rowは `T_i=[[A_i,0],[r_i,0]]` である。
先行blowdownのuniversal-quotient同型からH0(N*⊗A)=H0(T_P2⊗I_p)に一致し、追加cofactor仮定ではない。
net base b=(a,b,c)でrowsを(U_i,V_i)=(a,b)A_i+c r_iとする。

dual determinantを次の別chartで再導出した。
`y=(1,x,z), q=(−v*x−z,v,1), f_i=q*T_i*y, f=f1/f2`。
Zへのprojection bracketをFとして、kerdfから
`{x,v}=−F f_z/f_v, {z,v}=F f_x/f_v`。
base `(a,b)=(−v*x−z,v)` のbracketは `F(v f_z−f_x)/f_v`。
`f_v` numeratorはdet(y,T1y,T2y)=±F。
`D=∂x−v∂z`はqを固定し `y×Dy=q` だから、残Wronskianはdet(q,qT1,qT2)となる。
従って全identityは `tau=c Omega`、`Omega=det((U1,V1),(U2,V2))`。
これは一係数点の検査でなく、全matricesの有理恒等式であり、polynomialに戻る。

conic上のtwo rootsからD'lineの全quadratic coefficientsは
nodeで `(V1V2,−U1V2−U2V1,U1U2)`、cuspで `(V0²,−2U0V0,U0²)`。
primitive original rowsからこれらのgcdは一。
degreeが全て二なので、余分なbase rational multiplierはconstantしかない。
nodeの両rowsとcuspのdenominatorだけにprimitivityを課し、cuspの第二rowへ無断に課していない。

## 4. flat normalization、c|Omega の全補集合、L恒等式

Poisson submersionのsymplectic-base開におけるHamiltonian horizontal distributionはflat。
generic fiber上のregular vector fieldsはconicを保ち、Aut(P2,conic)→Aut(conic)のfaithful性から
conicのratio fを固定するconnectionで一意となる。
Sym²L frameでtrivialized vertical Poissonの係数kはJacobiによりnonzero constant。
conicはdet(L)²、projective fiber volumeはdet(L)³で変わるため、元のvertical cubicは
`(k/Omega)Q D'_fiber` である。

c|Omegaを除外する前にOmega_L≠0を仮定していない点を確認した。
Omega=c^m Omegatilde,m≥1なら、base inverse two-formはcにm+1位pole。
actual dense Y chartではtrivial fiber conic Q=c*v、flat conic local equationはOmega² Qのunit倍である。
kappa≠0によりcollapsed pointはstandard D'lineの外なので、flat fiber formはそこで
`dQ0/Q0 wedge regular-one-form + regular`。
flat projective coordinatesはunitを分母としてregularにpullbackするので、cに高々log単極だけ。
m+1≥2のbase poleと相殺できず、全c|Omega sourceが空となる。
任意m、全zero/special coefficientsを含むdivisorごとの矛盾である。

c∤Omegaを得た後、p^*Lはactual E又はD'のdivisor componentでないのでinverseはgenericそこでregular。
base residue db/Omega_Lとfiber residueを相殺し、全homogeneous polynomial identity

`Omega_L W(k1,k2)=k k1k2` （node）、
`Omega_L W(k1,k0)=k k0²` （cusp）

を得る。root∞、有限kappa zero、重根を除いたままにせず、全Lへpolynomial延長する。
generic parameterからspecial parameterへの外挿ではない。

## 5. Nodeの全binary分類とtrace係数の再導出

比例k_iはW0、kappa≠0に反する。
共通simple rootの場合にW divisibilityだけでは不十分である点を監査した。
**両k_i(p)=0ならrow(p)Aiはともにrow(p)と平行、従ってOmega_L(p)=0** を必ず使用する。
共通simple/simpleではordW≥2、ordOmega≥1となり右辺ord2に反する。
共通root無しなら各rootにordOmega=1が必要なので、degree二から両k_iは相異double roots。
残りはdouble/simple shared rootの31だけである。

tracesを含む二正常形を別に計算した。
22では `k1=s²,k2=t²` として
`A1=[[lambda1,1],[0,lambda1]], A2=[[lambda2,0],[−1,lambda2]]` を取ると
`Omega_L=−lambda2*s²+s*t+lambda1*t²`。
binary identityはlambda1=lambda2=0を強制し、両Aiは純nonzero nilpotent rank一。
31では `k1=s²,k2=s*t` として
`A1=[[lambda1,1],[0,lambda1]], A2=diag(lambda2,lambda2+1)`。
`Omega_L=s(lambda1*t−lambda2*s)`、W=2s²からlambda2=0、lambda1≠0。
全coordinate/shear形ではproducerのinvertible Jordan/rank一semisimpleと一致する。
root amplitudesはnonzero scalesで戻せるため、これはsource sampleでない全正常形である。

両型でOmega_Lは二つの相異単根。従ってglobal Omegaもreducedかつc因子無し。
ここで初めてOmegaのsimple residueを使う。

## 6. inverse residueと、producer §8 の一文の限定

c≠0でactual inverseはgeneric Omega×P2にregular。
flat conic equationのpullbackはOmega² Qなのでfiber formのresidueはnonzero constant*dlog f0、
f0はrank一L|Omegaのrow ratio。
base inverse residue eta_Bとの相殺は `eta_B=nonzero_constant*dlog f0` を与える。

**scope補正。** producer §8 の「eta_BのpolesはL∩Omegaだけ」は、Omegaがreducibleな段階では一般に正しくない。
Omega=ell1 ell2なら、各component上のresidueはell1∩ell2にもpoleを持ち得る。
著者もこの書き過ぎを確認し、producer bytesを変更せず本receiptで限定することに同意した。
正しい使い方は、reducible段階ではeta_Bがgeneric非零であることだけを使い、
**smooth irreducible conicを得た後に限って** poles⊂L∩Omegaを使用する、である。
以下の結論・全cover・dimensionには変更を要しない。

Omegaが二linesならmatrix image rank二以下では各componentのrank一imageが一点でf0constant。
rank三ならplaneとSegre quadricの交わりはopposite rulings二本で、一方でf0constant。
これはeta_B generic非零に反する。components交点のpole位置をこの排除に使っていない。
従ってOmegaはsmooth conic。matrix-plane embeddingとSegreの(1,1) curveによりf0はdegree一。
そのzero/poleは二rowsのzero pointsで、互いに異なる。
この段階ではeta_BのpolesはL∩Omegaだけなので、両row-zero点はL上にある。
31型のinvertible AiはL上row-zeroを持たず、全Ui,Viの一zeroはc≠0となって矛盾。
従って31の全sourceは空。

## 7. Cusp、真のfibers、全上界

cuspではOmega reducedを使用しない。
k0にsimple rootが二つあれば、common rootではdetrowsの追加zeroに矛盾、noncommonでは
各々Omega order二が必要でdegree二に反する。
従ってk0はdouble root。k1がその点で非零ならOmega order三が必要、double共通なら比例するので不可能。
従ってk1はその点にsimple rootを持つ。
全tracesを残したmatrix係数からdenominator A0は純nonzero nilpotent rank一、
numerator A1はそのfixed pointを保つ全三次元linear spaceの中に入る。
第二rowの非primitive境界、Omegaのdouble line、全lower rowsを保持したdimension上界である。

nodeはAi nilcone各二＋r_i各二の八parameter。
真のroot scales二はprojective pi_Eを保ち、global restriction injectivityによりprojective piを保つ。
free actionなので固定source像≤6。
cuspはA0二＋A1三＋r0,r1四で九parameter。
denominator/numerator scalesとnumeratorへのdenominator加法の三次元free actionが実projective pi fiberに入り、像≤6。
nodeに任意GL2 mixingを商として引いていない。
Aut(S)六はPGL5第一因子でfreeな対角作用を与え、全像≤24+6−6=24。

この全source proofを上のscope補正付きで受理する。
先行四notesの独立性の分担は§1の通りであり、新数値検証や全分類の完了を偽って主張しない。
