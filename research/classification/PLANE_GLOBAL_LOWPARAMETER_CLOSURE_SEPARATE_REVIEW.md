# rank零plane・globalVの全二十一低parameter源：独立閉鎖証明

2026-11-07、主義務G。producerと別作者による全source幾何及びFraction算術の独立確認。**PASS**。

`PLANE_RANK_ZERO_GLOBAL_VECTOR_LOW_PARAMETER_CLOSURE.md` の二十一whole源について、
全higher Poisson tensorのambient projective像次元が24以下であることを受理する。
ここでwholeとは二つのfirst-jet noteで**指定した各全source**を指す。
generic二・三parameter源、旧m二外の全復活source、全plane、全G/Eを尽くしたという意味ではない。
normalnil b非零の全mへの補足は別の
`PLANE_NILPOTENT_AFFINE_ALL_MULTIPLICITIES_SEPARATE_REVIEW.md` で受理済みである。

## 1. immutable入力と有限job

全packet入口は `inputs/plane_rank_zero_global_vector_lowparameter_packet_manifest_20261107.json`。
二十三件のうちmutable runner一件とjobs八件は、現在のfileではなく
`inputs/plane_rank_zero_global_lowparameter_frozen_inputs_20261107.json` の保存bytes/hashから読んだ。
その九snapshotと残る十四immutable fileの全hashを独立に確認した。
後続two-parameter modeを足した現runner又は現jobsへ旧証明書をbindし直していない。

独立checkerは `code/classification_plane_global_lowparameter_closure_review.py`、結果は
`reports/CLASSIFICATION_PLANE_GLOBAL_LOWPARAMETER_CLOSURE_REVIEW_20261107.json`。
全入力・libraryのhashと240秒/1536MiBの有限契約は
`jobs/plane_global_lowparameter_closure_review.json` に固定した。
gate PASS後、13.233秒、peak RSS49,252KiBで全検査PASS。
初回は先行libraryのpolynomial parserが `tau**(-1)` というrational kernelの表記を受け付けず停止した。
新reviewer内の専用Fraction(τ) parserで負の整数冪を明示的に扱い、契約を再gateして検証した。
これは数学上の反例又はproducerの変更ではない。

producer、producer helper、SymPy、CAS又はnative GBをimport・実行していない。
先行の別作者独立Fraction・疎多項式・Schouten算術だけを再利用し、全行列をtensorから別に再構成した。
旧LOWPARAMの `all_passed=false`、compatibilityでのT-square印字失敗及びnative探索記録は不変である。
今回の最終証明はnative ideal dimension、GB完了、ideal lifts又は初期idealの等式に依存しない。

## 2. first symbolの全幾何と正しいbase次元

normal半単純型では、trace-free base linear fields T,Vを用いて

`Pi0=T∧En+V∧(z∂z−w∂w)+kzw∂z∧∂w`, `[T,V]=0`。

m二はdet(Ebase,T,V)非零と同値で、span(T,V)は独立な二次元abelian代数である。
regular matrixを含めばそのcentralizer自身となり、Cartan、mixed〈H,N〉、regular-nil〈J,J²〉の三型。
含まない場合、nonzero semisimple要素があれば別のcommuting nonscalar要素とのgeneric和がregularとなって矛盾する。
従って全要素がrank一nilpotentであり、rank一matrixの和のrank条件からcommon image又はcommon kernelの二型だけ。
common imageはpointwise proportionalでm二に反し、common kernel〈E02,E12〉が残る。
Cartanは既済の全MPS定理、残る三型は公開noteのactual normalizerと非零tensor scale及び逆で列挙される。

normalnil b非零の全source及び共鳴は先行独立レビューで確認済みで、今回も全operatorを再構成した。
b零では `[T,V]=0,K=az²` であり、Vを冪零に限定しない。
actual normal diagonal changeとtensor scalarはT,Vを独立にscale出来るので、
二つのfirst-jet noteのmixed、regular-nil、common-kernel全normal formを使える。
全零係数、square-line、semisimple-line、resonant rowを別に保ち、generic二parameter源を一つの標本へ縮めていない。

full affine symbol base BはsubalgebraのPGL3 orbit、ordered pair四、normal行列line、Kから数える。
半単純normalではmixed orbit五、nil orbit四、common orbit二にnormal line二、K一を加え、
一般sourceは順に12,11,9。Vがsquare lineなら一減、K零なら更に一減する。
冪零normal b零ではnormal line一なので一般mixed11、regular-nil10、common8。
T又はVが固有のsemisimple/square lineへ入れば各一減、K零も一減する。
b非零のB11、固定ρ又は追加共鳴orbitのB10は先行独立全源証明による。
これらは消費した座標正規化を戻した全symbol次元であり、normal-form parameter数ではない。

## 3. 全二十一tensor、全350×30/9と全wall

checkerはH=diag(1,1,−2)、N=E01、J=E01+E12、J²=E02から各T,Vを別に作り、
各sourceの全150係数、正しいEuler-contracted実Pfaffianを再構成した。
全30/9 symbolic source及びraw36=canonical30+Euler6のrankを独立に確認した。
各Pi0から全350行のD1,D2、六つのactual linear shearの全30×6行列を再構成し、保存全entryを照合した。
各nonzero maximal minorとcomplete kernelの全行代入零・独立性を確認した。
全一parameter源の選択minor積の完全一次分解と、全rational wallでの別minor/kernelも認証した。
未確認の代数的wallは無い。

全higher解ではPi1∈ker D1、固定Pi1上のPi2は同一D2のinhomogeneous fiberなので
`d≤nullity D1+nullity D2`。plane六とtensor射影化−一を戻した像はB+d+5以下。
fixed rootではBを一減する。
この全算術で既に閉じる十六whole modelsの上界は、元noteの表と同じく
20,22,22,21,21,21,20,24,20,23,21,22,21,21,20,19である。
全generic openと全wallを検査したので、generic rankだけの外挿ではない。

未閉だった二つの一parameter源のK非零部分も全21以下。
K零では両者ともB9、D1 kernel八、D2 kernel三で、線形和だけなら25となる。
最後の三fixed sourceの線形和は28,27,27であり、それをそのままG boundと呼ばない。

## 4. 二つのK零源を全cokernel一cutで閉じる

normalss T=J,V=J²,K0、normalnil T=J+J²,V=J,K0の全D1 kernel八方向をu0,…,u7で保つ。
その全Pi1からQ2=Jac(Pi1)のprojective全350係数を直接再構成し、保存全式に一致した。
全D2を消す次のfunctionalを全九columnで検証した。

* ss: `(2/3)e51−e61+e132`、Q2への値は `−4(u2−u7)²/9`。
* nil: `(1/3)e51−(2/9)e52+(2/9)e53−(1/3)e61+(1/9)e62−(2/27)e63+e131`。
  Q2への値は `−2(9u1²−15u1u3+36u1u5−12u1u7+5u3²−30u3u5+5u3u7+36u5²−24u5u7−u7²)/9`。

両方非零polynomialなので、全Pi1必要源のdimensionは七以下。
D2のnonzero六minorと逆も再認証し、全Pi2 fiberは三以下。
従って全像は `9+7+3+6−1=24` 以下である。
特殊uを選んだ空性ではなく、全八係数への必要な一条件である。

## 5. actual有限shearと全weight二の解空間

残るss common K0、nil common K0、nil T-square K0では、
actual linear change x→x+C(z,w)の有限作用のfirst-weightは正確にPi1→Pi1+S yとなる。
normal weight一のgeneratorなので、first-weightに二次以上の作用は無く、全higher差はweight二の全九係数へ入る。
逆はminusである。

独立なS列と保存Rを合わせた全kernel basisについて、nonzero basis minorと全inverseを検証した。
全D1がRを消し、そのrankと列数がcomplete kernel次元に一致するため、
`ker D1=im S ⊕ im R` を全source上で得る。
ss/nil commonはshear rank四とslice八、T-squareはrank五とslice六。
このactual有限群から、全higher sourceはslice sourceにshear rankを加えた値以下となる。

Rの全Pi1からQ2全350を独立に作り、D2のconstant六minorで六係数を解いた。
**complete kernel D2の三係数v0,v1,v2を全て戻した**。
各cokernel functionalが全D2を消し、Q2への値が保存generatorに一致することを検証した。
さらに全350残差と表示generatorsのQQ-linear rowspacesが一致した。
全式はquadraticなので、これは必要な選択式だけでなく全weight二compatibility idealの等式を与える。

## 6. common二型の全linear branchesとheight

ss commonの全compatibilityは非零定数を除き

`u7(2u2+u3+2u5)=0`, `u4(−2u0+u1−2u6)=0`。

従って全reduced台は、u7又は第一linear form、u4又は第二linear formを零とする四linear spacesの和である。
全零境界と全intersectionを含む。各branchはu六と全v三の九変数。
その全weight二解を戻し、元の全350 weight三及び全350 weight四を独立展開して全entryで照合した。

各branchのrow65=F、row66=GのgcdはC上でも1である。
この独立checkerではproducerの多変数gcd結果を再利用せず、次のexact certificateを用いた。
Fはu0に関して三次で、そのleading coefficientは非零**定数**である。
従ってC[他変数][u0]の任意の非定数共通因子は、unitを除きu0についてmonicで正次数となる。
Fのcontentがunitであるためu0に依らない共通因子も無い。
他変数をどの定数へspecializeしても、そのmonic因子の次数は落ちず、特殊化F,G双方を割る。
一方、u1,u2,u6を1、他のu0以外を0としたQQ上の二つの一変数polynomialについて、
独立Euclidean division chainのgcdは全四branchで非零定数だった。
その全polynomial、quotient、remainderを新reportへ保存した。
従って全C-sourceに共通因子は存在しない。
これは任意の点評価からgeneric rankを推論する操作ではなく、monic因子の次数保存で完結する全体証明である。

polynomial UFDのheight一primeはprincipalだから、gcd1の二式を含む全Jacobi idealはheight二以上。
unit ideal又は空とは主張しない。
各sliceはdimension七以下で、全像は `7+shear4+symbol8+plane6−1=24` 以下。

nil commonの全compatibilityは

`u7(u3+u7)=0`, `u7(−u1+2u2+3u5)−u3u5=0`。

u7零ではu3又はu5が零、u7非零ではu3=−u7,u1=2u2+4u5となる。
従って全台は、この三linear spacesの閉和で尽くされ、u7零への閉境界も含む。
各branchのu六+v三について同じ全700係数を再構成し、weight三row65は全て非零polynomialだった。
よってsliceはdimension八以下、全像は `8+4+7+6−1=24` 以下となる。

## 7. T-squareの全reduced supportとoriginal一cut

全七compatibility generatorsをg0,…,g6とし、L=−u2+3u3+u4、d=u2−u4とする。
独立に得た全polynomialsについて、次のliteral identitiesを検証した。

`g1=−12u5²`, `g3=u3L/3`, `g5=14u0u5/3+L(2d+9u3)/6`,

`L²=(15/2)u3L−(1/2)L(2d+9u3)`。

全reduced点でu5=0、u3L=0、L(2d+9u3)=0だからL=0。
逆にu5=L=0を全七generatorへ代入すると全て恒等零であることも検査した。
従って全reduced supportはu5=0,u2=3u3+u4という単一の四次元linear spaceである。
scheme-theoretic radical又はnilpotence指数を追加仮定していない。

この全spaceでP2の全inhomogeneous解とkernel三係数を戻し、全weight三・四700係数を再構成した。
元のweight三row65は正確に `u3(2u0²−v1)` であり、七変数上の非零polynomial。
従ってslice sourceは六次元以下、全像は `6+5+8+6−1=24` 以下である。

## 8. 最終射程

以上で二十一whole低parameter源の全generic open、全wall、五つのK零frontierを閉じた。
旧十六good源の証明まで未解決扱いに戻さず、同時に二十一を越えるgeneric源へ拡張しない。
normalss generic MIX三parameter、MIXH/MIXN/NILREG二parameter、normalnil b零Cartan三parameterとgeneric mixed二parameterは別義務。
全mでの共通入口は成立しても、その全復活normal formsをこのpacketだけで尽くさない。
新二parameter四源・hyperbola packet、成分性・全E帰属、全G及び主分類は本reviewの範囲外である。
909主原稿/PDF、global台帳、Gitは変更していない。
