# Q_b自己旗stratum上のactual CA/CB交叉：nil module Z4の全環

2026-11-03。先行P_aのfull-D型をラベル交換せず、非対称な実QQ(b) embeddingと全ambient cubic gateを使用する。
結論はQ1および非空unit openとそのG-saturation上の形式環についてであり、b=0のB18原点や全bad locusを分類しない。

## 既に受理された配置と、core idealの正当化

CAはsmooth29、Nはsmooth28。二つのCB branchをB1,B2とするとB2∩CA=Nはcleanである。
B1∩CAの全相対完成環はZ4=P[u,v]/(u²,v²)、P=C[[s20]]。
二CB branchのmeeting C2=P[v]/v²への実quotientはuを殺す形で、kernelがisotropicであることもQQ(b)全source計算で確認済み。
従ってCAのregular29環BでN=(f)とし、座標を整えるとactual Z4 idealは

    J_Z=(z1,...,z7,f²,rho²),  B=P[[z7,f,rho]]。

二CB branchの全fiber-productには、共通rho、branch方向u8/v8、tauを使う全81関係がある。
ここではbranch1のidealを(v8,tau)、branch2のidealを(u8,tau−rho²)とする。
Nはbranch2に入り、CAのN横断方向はbranch1内でbranch2を外れるので、一つのu_iのCA制限はf×unit。
全u_iとtau−rho²のCA制限はfの倍数である。
従って81関係のCA制限は正確に

    f*( v1|CA,...,v8|CA,tau|CA )

を生成する。

この括弧を全J_Zと同定するgateを省略しない。
actual CBの残りcommon ambient equationsの二normal quadratic blockは、全38 tangentへの厳密projectionで零である。
よってZ4の七つのlinear unitと二つのregular quadricsは上記coreの制限から得られる。
degree3 spanの単元性とNakayamaにより、このcore ideal自体が全J_Zに等しい。
従ってactual D ideal Qには**全idealとして**Q0=f J_Zが入る。
単なるfiber productの非flat basechange exactnessを仮定したのではない。

## extra ambient equationsを全て除く最後の一bit

D_red=Nであり、またDはZ4を含むため

    f J_Z ⊂ Q ⊂ (f)∩J_Z。

右と左の差はP上でf²,f²rhoの二classを持つmoduleである。
全S上のactual D quadratic-zeroはf² classを除き、残りはfree rank1 module P*f²rhoのsubmoduleに入る。
ここでextra ambient equationsを無視することはできない。

同じ実29 target baseを使う三sourceは全てétaleであり、CBの二sheet環はE⊕tI（t²=1）となる。
**全121 ambient odd coordinate**のb-vector、CA sectionのc-vectorを保持したexact colon recipeを使用する。
Iは8linear+1quadraticのregular sequenceで生成される。
全b-vectorのregular-sequence divisionを実三sourceのthird-order inverseで計算し、全150行のsyzygy identityをQQ(b)上で確認した。

8linear generatorのcofactorはdegree2まで必要なので、実b sourceはdegree3まで保持した。
quadratic generatorのcofactorをdegree1までで止めてよい理由は、対応するactual c_quadのconstant termが厳密に零だからである。
これもcertificate publisherで明示的に確認する。
したがって未計算のdegree2 cofactorやfourth-order bを落として目的のdegree2 defectを変えることはない。

全150のdefect quadraticはactual Z4 idealへ入り、f²rhoの係数は**QQ(b)で全て恒等的に零**だった。
GL×bは実S20 stratumをdominantに覆うので、その全generic moduleも零である。
Pはdomainであり、残差はfree rank1 P-moduleのsubmoduleだから、genericに零のtorsionを残す余地もない。
全higher normal powersは既にJ_Zへ入るため、blind fourth-jet計算は不要。
ゆえに

    Q = f J_Z = (f z1,...,f z7,f³,f rho²)

が全completed idealとして証明された。

## 全ring、nilpotence、primary structure

従ってactual Dの全形式環は

    C[[s1,...,s20,z1,...,z7,f,rho]] /
       (f z1,...,f z7,f³,f rho²)。

nilradicalは(f)で、そのB-moduleは

    (f) ≅ B/J_Z ≅ P[f,rho]/(f²,rho²)

となる。P上free rank4、nilradical²はP[rho]/rho²に同型なrank2、nilradical³=0。
これはP_a側のsquare-zero rank1 nil moduleとは異なる。

irredundant primary decompositionは

    (f z7,f³,f rho²) = (f) ∩ (z7,f³,rho²)。

第二idealはprime (z7,f,rho)にprimaryで、generic primary quotientはlength6（f³,rho²）。
**nil moduleのlength4とembedded primary quotientのlength6を混同しない。**
associated primesはこの二つだけ。
socle要素f²rhoは非零で全九normal変数に殺されるので、normal coreのdepthは0。
Pの20smooth variablesを戻すとcone dim28/edim29/depth20/multiplicity1。
射影ではdim27/edim28/depth19/multiplicity1、embedded supportはcone20/projective19である。
reduced support Nはここで正則であり、full DはS20上で非被約・非CMである。

## global embedded primeへの帰結と範囲

S_Qは実GL×b像の既約閉包であり、上のunit openにおけるembedded associated primeはそのalgebraic closureへdescendする。
忠実flat completionでのassociated-prime contractionと実G-saturationを用いる。
先行S_Pではwhole CA flag count2、CB count1であり、本S_Qでは逆にCA1、CB2。
両closureが同じ既約stratumなら、その二つの非空generic openが交わってflag countが矛盾する。
従ってactual Dには**二つの相異なる既知embedded prime closure**がある。
これは全embedded-prime listのexhaustionや新しい全Poisson componentの追加ではない。
原点B18の相対CA/CB配置、全D型、両stratumの他の境界点は別課題として残す。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CB_self_pair_actual_D_full_ring_20261103_v63.py --output <new-output-directory>

全source計算は別のcubic gate packetへ保存されている。
このproducerは数学的applicationの依存を固定し、全monomial primary intersection・nil basis・socleを独立に確認する。
