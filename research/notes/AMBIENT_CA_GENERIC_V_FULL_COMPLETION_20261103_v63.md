# CAのgeneric V境界：actual成分の全完成環

2026-11-03。本稿は `AMBIENT_CA_V1_ODD_ORDER_RESULT_20261103_v63.md` の適用鎖を独立監査PASSの後に適用するcompanionである。
再現プログラムはその監査reportを必須引数とし、PASSと全依存ハッシュを確認するまで結果を出さない。
V1とその稠密G-orbitを扱う。B18原点のCA full ideal、CB、actual CA/CB scheme intersectionの型は主張しない。

## 物理的sourceから全成分への適用

CAのproper flag modelは既知成分の全proper像であり、非零点のsourceは滑らかな29次元である。
P0における全proper flag fiberは一台長2と既に証明されている。
従ってP0を含むG不変開集合でproper写像はfiniteで、全fiber長は2以下である。

全critical schemeのmonic V-axisはP0を含む代数曲線である。
そのG飽和閉包Zは既知のcritical schemeに含まれ、P0におけるその次元は20である。
V1のGL orbit微分rankは20なので、Zは20次元、V1 orbitはZの稠密開軌道となる。
上のfinite開集合はZに会い、G不変性によりV1 orbit全体を含む。
V1 source微分はcorank1なのでfiber長は2以上となり、全fiberは一台長2である。
ここまでは滑らかなsourceからのfinite写像を得たものであり、proper性だけからbirationalityを仮定しない。

V1の全150係数によるSchur critical cotangentはrank9である。
一方critical schemeは20次元orbitを含むので、その局所環は正則20次元となる。
実28 target座標と、余剰一target座標を使うquadratic base Rは正則29次元であり、

    B = R[t]/(t²−nu),       A = R ⊕ I*t ⊂ B

がfinite sourceと実成分の完成像である。
critical quotientは B/(t,I)=R/(nu,I) だから、J=(nu,I) はRの正則9生成idealである。

八つのactual regular odd constraintを解いてから、残る出力の三次係数
800/27 と −40/27 が非零であることを厳密QQで確認した。
さらに全150 actual odd coefficientのfirstjetはrank9で、Jの全九classを張る。
よって J=I+m_R J であり、有限R-module J/IへのNakayamaにより **I=J** である。
Rのformal座標を選べば

    R = C[[g20,nu,s1,...,s8]],       I=(nu,s1,...,s8)

と書ける。Iは非零なので、非零c∈Iを取ればt=(c*t)/cがFrac(A)に入る。
従ってfinite normal ring BはAと同じ分数体を持ち、ここで初めてBが実normalizationであると結論する。
この順序はglobal proper flag modelだけから未証明のbirationalityを持ち込まない。
これはraw Y固定sliceのcuspからの類推ではない。
G20 orbit方向をformal sliceから戻す操作は成分のG不変性を用いる。
チェックポイントのGL列10は全150係数で点P自身と一致するので、Poisson scalar方向もこのorbitに含まれる。
非零scalarを除くprojective局所環ではg変数を19個に減らす。

## 全完成環と完全81関係

従ってV1でのcone完成環は

    A = C[[g20,s8,t²,s1*t,...,s8*t,t³]]
      ⊂ B=C[[g20,s8,t]].

core baseを b=(nu,s1,...,s8)、九odd targetをY_i=b_i*tと取る。
全target kernelは次の81関係で生成される。

* 36個のKoszul関係 b_j Y_i−b_i Y_j、i<j。
* 45個の積関係 Y_i Y_j−nu b_i b_j、i≤j。

積関係で任意の式はR+ΣR Y_iに還元される。
その像が零なら、BのR-basis (1,t) により定数項とΣc_i b_iが零となる。
bは正則列なのでその第一syzygyはKoszul関係で全て生成される。
従ってpullback零の検査だけでなく、これが全kernelである。
再現プログラムは別のMacaulay2 elimination用scriptも保存するが、本checkpointでは実行しない。
全kernelの証明は上の正則列Koszul argumentによるものであり、未実行のGBを証拠に数えない。

## 不変量と正確な局所bad locus

coneではdim29、embedding dimension38、depth21、Hilbert–Samuel multiplicity2である。
projective版はdim28、embedding dimension37、depth20、multiplicity2である。
depthはA=R⊕Iと正則列九生成idealのKoszul resolutionから得られる。
odd generatorの二乗はm_R³に入るのでm_R Aはm_Aのreductionであり、multiplicityはR-rank2である。

normalizer conductorは **(t²,s1,...,s8)B**、target conductorは
**(t²,s1,...,s8,t³,s1*t,...,s8*t)A** である。
その商A/conductorは正則C[[g20]]、B/AはC[[g20]]*tである。
Aはdomainかつanalytically irreducibleであり、embedded associated primeを持たない。
このformal germではsingular、nonnormal、non-CM locusは全て、この滑らかなconductor supportに等しい。
conductorの外ではA=Bであり、support上では正則九生成idealIが非principalであるためAはCMではない。

t²,t³はAに入り、tはAに入らない。
従ってB=A[t]はelementary subintegralなnormalizationで、finite universal homeomorphismとなる。
seminormalizationはBであり、Aはseminormalではない。
これはreduced base上の異なるnormalization二点を排除するが、長2のscheme fiberをreducedとする主張ではない。

この結論をB18原点へ延ばして、全cubicをU³と同定してはいけない。
U-axisの三次witnessでは、nu*Vを経由するUV²のような混合項は検出できない。
CBのgeneric V計算と、原点でのCA/CB full odd idealは未解決のままである。
