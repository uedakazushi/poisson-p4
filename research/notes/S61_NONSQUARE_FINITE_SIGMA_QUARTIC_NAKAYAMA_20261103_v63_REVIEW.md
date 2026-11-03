# 有限Sigma上quad78の独立レビュー

2026-11-03。author packetの結論 **元raw quadratic-cofactor mapは有限closed Sigma全域でrank78** を、別のFraction係数演算・逆monomial順のnative matrix構築と一回の独自rankで点検する。全開quad78、実C122 rank64、fullquartic全射は受理しない。

## 独自算術検査

元の804quad係数を、元六light chartからoriginald→P/a・a-clearing shifts `[0,2,2,1]` により独自再構成する。Q2primeへのrow変換は使わない。保存21×4 right inverseについて全16のQQ恒等式 `QL=gI4, g=-a²P²J` を独自Fraction dictionary arithmeticで検証し、全2461cofactor係数のp=32003分母がunitであることを確認する。

六Koszul row matrix K（6×84）のright inverseは、列(a,b),a<bのblock aへL[:,b]/gだけを置くordered choiceである。独自に求めたQLから全36のKR成分がI6となることを照合する。K*M=0は、同じquadratic-cofactor indexにおけるq_j*q_i−q_i*q_jの全係数が零という普遍多項式恒等式である。native側の独自monomial構築がそのindexを固定する。

既存の独立finite-Sigma reviewで証明された28次元special algebraのP,e operatorsをhashで固定して再利用する。新nativeは28×28operator arithmeticもFLINTで行い、authorとは逆のdegree2/4 monomial順で実2352×3528 blockを構築する。元sourceから作った4×21係数だけを入力し、authorのraw24行抽出codeを呼ばない。

authorが保存したblockとの全8,297,856成分を、反転monomial permutationを通して比較する。この比較は座標順・operator orientationも含む。その後、一回だけ独自FLINT rankを計算し2184を確認する。大きなinverse、別field、QQ GBは作らない。

## Orientationと有限Nakayama

標準multiplication operatorはoutput-row/input-columnなので、blockは **M:S^126→S^84** を表す。六KoszulはK:S^84→S^6であり、K*M=0。逆向きのgenerator mapをこのblockと同一視していない。

既存proper finite model SはV=Z_(32003)上finiteで、FULLUはS全体のunit。従ってgもunit、Q,Lの数値分母もp-unitなので、上のKR=I6はS上で定義される。H=kerKはS^84のdirect summand/projective rank78であり、V上finite。そのsplitのためH/pH=ker(Kbar)となる。Sのflatnessを仮定しない。

special algebra A=S/pSのFp次元は28。従ってdim Hbar=78*28=2184。Mbarの独立検証rankが同じで、imageMbar⊂HbarなのでMbarはHbarへ全射。C=H/imageMはfinite V-moduleでC/pC=0。NakayamaによりC=0である。

したがってS上でMはHへ全射、cokerMはS^6と同型。Hはprojectiveなのでsurjectionはsplitし、任意のQQ/complex fiberでrankM=78となる。minor idealやArtin socleを介した推論も不要である。nonproper openのfinite-field rankを直接QQへ移していない。

全Sigmaでの元ell chartと既証ZW19を合わせるdegree3-colon argumentにより、original light146・weighted152もこの有限Sigma上で従う。これは他の有限例外や全coreを含む主張ではない。

## 保存と再現

code `s61_nonsquare_finite_sigma_quartic_review_20261103_v63.py/.cpp` はauthor implementationをimportしない。既に独立受理されたfinite modelのoperatorsだけを固定premiseとして再利用し、元quad・QQ split・新blockを別に作る。生成ELFを保存せず、33MBblockは元のignored workspace sourceをhashで固定する。新しいcompact4quad入力とnative receiptを保存する。

再現は既存Cloud activation後 `python code/s61_nonsquare_finite_sigma_quartic_review_20261103_v63.py --output <new-output>`。一回のrankを含め60秒・2GiB上限とする。
