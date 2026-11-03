# 指定 Q の実交叉：全有限表示・主平方零 nilradical・11次元埋込み素イデアル

2026年11月3日。対象は S61 の二つの既知の **被約成分** CA,CB の実際の scheme-theoretic
intersection D=CA∩CB である。指定点は継続している v24 Q(a=1)であり、a=2/9の別点ではない。
全P4におけるgeneric Poisson rank≤2だけに支えられた全域境界はこの定理の対象外である。

## 定理

指定 Q の近傍で次が成り立つ。

1. D の全完成局所環は、下記49変数・218個の明示的な多項式関係の商に、
   七つの独立な smooth orbit variables を加えたものである。
   正規化同士の tensor product ではなく、実際の二つの target image ideals の和である。
2. nilradical は主 ideal (F)で、F≠0、F²=0である。従って正確な冪零次数は2である。
   被約商は既存 N の完成局所環 A7 と一致する。
3. nilradical の全加群は

       (F)≅C[[u1,…,u11]]*F

   で、その零化 ideal は滑らかな11次元部分の prime ideal である。
   これは共通の10次元非正規層より一方向大きい。
4. D の associated primes は正確に二つである。
   N の minimal prime と、この11次元部分の embedded prime である。
5. S61 の完成局所環として dim D=18、depth D=11である。
   値・perspective の九つの smooth directions を戻すと dim27、depth20となる。

環の直和分裂 `D≅A7⋉C[[u11]]F` は主張しない。
全環は保存した実 target ideal の有限表示で決まり、nilradical とその全作用は以下で決まる。

## 1. 同じ target K を使って全 normalizer fiber を多項式で戻す

既存の臨界一次ジェット slice

    τ=p,A=0,B=C=1,D=−p/2,F_base=s,  旗=(0:1:0:p)

を使う。ここでF_baseは本ノートのnilpotent Fとは別の基底係数である。
raw CA/CB の全61出力から、共通の十個の実q座標をK=(v0,…,v9)に選ぶ。
その行は45個のq内の0始まりindexで

    (0,19,21,24,25,27,28,32,33,34)

である。両族の source fiber十係数からKへの全係数行列Lγ(p,s)は

    det L_CA=det L_CB=−72

を恒等的に満たす。従って十のsource fiber係数はK,p,sの **多項式** として全て戻せる。
近傍でのみ分母を除いた近似ではない。

残る35個のqをz0,…,z34とし、a=p²,b=p³,c=spと書く。
共通の実target ambientを

    R=QQ[z0,…,z34,b,c,v0,…,v9,s,a]

という49変数環とする。各CA/CBの実q出力をpについて偶奇に分けると正確に

    z=Eγ(K,s,a)+p*Oγ(K,s,a)

となる。一次ジェットの残りの係数もa,b,c,sの明示多項式で戻るので、
この49変数表示は元の61座標の共通閉基底像方程式を除いた全表示である。
七次元Sp4因子の形式的product分解は既に受理された全slice定理を使う。

## 2. 各成分109式が全 kernel を生成する証明

Oγ(K,0,0)の係数行列は階数7である。
その七つの独立な線形形式をj2,…,j8とし、j0=a,j1=sと置く。
J=(j0,…,j8)はT=QQ[K,s,a]のregular sequence idealである。
残差をaとsで割る全多項式分解により、zの全35出力は

    z=Eγ+Cγ*(j2p,…,j8p)+b*Uγ+c*Vγ

と書ける。Cγは定数行列で階数7である。
従って28本の実target線形関係を除くと、残る成分像は正確に

    T[y0,…,y8]→T⊕Jp ⊂ T[p]/(p²−a),   yi↦ji*p

で与えられる。
この写像の **全 kernel** は

    ji*yj−jj*yi       (i<j),
    yi*yj−a*ji*jj     (i≤j)

で生成される。実際、二番目の関係によって任意の多項式はT+ΣTyiに戻る。
像で零なら定数部分は零であり、線形部分の係数はJのsyzygyである。
regular sequence idealの第一syzygiesはKoszul syzygiesで生成されるので、最初の関係で尽くされる。
従ってこれらは必要条件だけではなく全kernelである。

各族の式数は28+36+45=109である。
それぞれのtarget idealをIA,IBとすると、実交叉idealは

    ID=IA+IB

であり、保存した218式は実交叉全体の有限表示となる。
全式のsourceへのpullbackも恒等的に零であることをQQ上で確認した。

## 3. 一つの具体的な nilpotent と被約商

元の61座標を0始まりでu0,…,u60と書けば

    F=u0−4*u34+u51

と取れる。49変数環では

    F=z25+v0+2a.

CA sourceへのpullbackはrA、CB sourceへのpullbackは2rB−25p²/2である。
CA上のpFも実target多項式Gとして全式で構成した。

IAの像をF=0へ制限するとNとなり、そのtarget idealは

    IN=IA+(F,G)

である。ここでNの九個の実q座標とFをsource fiber座標に使う全Jacobian determinantも−72なので、
source fiber restriction rA=0 を全パラメータで正確に扱える。
さらにA_CA⊃R′,p²Bの全像記述から、rAで消えるtarget idealは(rA,rA*p)に厳密に生成される。
従って上のINはNの全prime idealである。

独立に実行したactual ideal計算は

    F mod ID≠0,
    G mod ID=0,
    F² mod ID=0,
    ID+(F)=IN

を与えた。最後の等式は全生成元を両方向へreduceして検査した。
INはprimeなので、これによりnilradicalは正確に(F)、平方零である。

## 4. 全 annihilator は滑らかな11次元 prime

同じactual target idealについてcolonを計算すると

    ID:(F)=Q_emb

となる。Q_embは49変数空間の45本の独立な線形形式が生成するprime idealであり、
残る自由座標は正確に

    v1,v4,v5,v7

である。colonと線形idealの一致を両方向の全ideal membershipで検査した。
従って七次元群因子を戻したnilradical加群は正確にC[[u11]]Fとなる。

p=s=0でこれらのfiber座標は

    v1=−k1,   v4=k0,   v5=2k2,   v7=3k5

である。
元の共通非正規層ではk1も零であり、自由係数はk0,k2,k5の三つであった。
しかし **k1*Fは非零**であり、k1²*Fも非零である。
従ってnilradicalを10次元非正規層だけに支えられた加群と置くことは誤りである。
埋込み層はk1方向を自由に含む11次元で、そのk1≠0部分は個々の成分が正規・滑らかな点も含む。

正確な完全列は

    0→R/Q_emb→R/ID→R/IN→0,

で、左の単射は1をFへ送る。
左項のassociated primeはQ_embだけ、右項はdomainでminimal prime INだけである。
一般のassociated-prime包含と、最小prime・実際の部分加群からの逆包含により

    Ass(R/ID)={IN,Q_emb}

となる。他のembedded primeはこの局所模型にはない。

完備化と七次元因子を戻した後、左項は正則11次元、右項A7はdepth11である。
local cohomology完全列でH^10(A7)=0、H^11(C[[u11]])≠0がDへ単射するので、depth D=11となる。

## 5. 計算記録と範囲

producerと全receiptは

* `code/ambient_intersection_exact_ideal_20261103_v63.py`
* `code/ambient_intersection_exact_ideal_finalize_20261103_v63.py`
* `reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json`
* `reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/report.json`
* 最終M2 script/log: `intersection_colon_and_reduction_v7.m2/.log`（同じreport directory）

にある。全109+109式、座標の全逆写像、G、colonの線形prime生成元、全SHAを保存した。
M2のProductOrder{35,14}ではactual IDのGBは実測0.04秒で完了した。
先行GRevLexは120秒上限で停止したため、変更したadapterを別名で保存した。

先行した二つの環分裂候補は全ideal equalityを満たさなかった。
それらを成功と扱わず、script/logと最終reportに保存した。
本定理には環分裂は不要で、実有限表示・nilradical・annihilator・associated primesが全て確定している。

S61とrank4開部分の範囲は `AMBIENT_INTERSECTION_SYMPLECTIC_SCOPE_20261103_v63.md` に従う。
この局所結果をgeneric Poisson rank≤2の全域境界や全Poisson分類へ広げない。
