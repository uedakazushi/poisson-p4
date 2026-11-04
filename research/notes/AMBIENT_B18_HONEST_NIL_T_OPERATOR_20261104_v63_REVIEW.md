# honest T_nil128×972：独立 source／形式函数／完全性監査

2026-11-04。元の零次数実環 D0 の nilideal N と reduced-normalizer cokernel C' を扱う。
S=k[[zA0,…,zA6,rA]] は同じ実 F9 inverse による ell=0 特殊化であり、rA は DVR parameter
ell ではない。新しい source値、jet、点、GB、823680 scalar mX 積は計算しない。

## 実sourceと全除算

受理済み radical complex は X:44×780、V:3×44、H1=N、H0=C' で完全である。
全640 Brad/P source と全140 SyzQrad sourceを保持する。targetはP^-1 V、source changeは
Z=diag(Res4(A^t),I8)。mu の中へ転置を持ち込まない。dh と dcl および四C-basis倍は
Qrad の文字通り列20〜23、28〜31で、Vのodd source符号とは混同しない。

独立な固定QQ arithmeticで、新gammaの四つの実商class
(gamma1,gamma2,gamma3,rA gamma3)が全C'barの基底となること、全zが零で、rAは最後二つの
間だけ作用することを確認した。従って全次数でsector idealsは(z7,r),(z7,r),(z7,r²)。
これは既受理の完全有限module作用からの結論であり、有限jetからの外挿ではない。
三Koszul differentialをwedge削除符号から再構築し、全列を照合する。

新mの132列成分と旧補正Hの72列成分について、全204 canonical decompositionを実函数の
添字まで確認する。各々の七prefix差はそのcoordinateで割り切れ、r-power部分は上記の
完全cokernelへの零像から割り切れる。旧CB-sectorではy=z+c r²の逆shiftを先に、正shiftを
最後に適用する。全1632 division slots（同一nodeの共有あり）を確認し、rやconductorの
逆函数を仮定しない。固定24残差は保存された八nil像の定数結合で全八coordinateが一致する。
これにより npre=Z nold H+L B は実函数として Vgamma npre=d1 を満たす。

## 単元証明と780 lift

Vgamma=d1 m より d1(m npre−I)=0。完全minimal Koszul resolutionのd2は全成分が極大ideal内
なので、M=m npreはM(0)=Iとなる。従って形式函数行列Mの逆は正則であり、n=npre M^-1は
mn=I、Vgamma n=d1 を満たす。scalar Mを評価したという主張ではない。

各sectorの元のmX列vについて、cycle identityをz=0に制限するとr^n v7(0,r)=0なので
v7(0,r)=0。ordered Hadamardでv7=ΣziHiと書き、wi=vi+r^nHiとする。
Σziwi=0に対し、pair(i,j)<7の係数を
(Ez+2)^-1(∂i wj−∂j wi)、pair(i,7)をHiと置く。δの規約はδ(ei∧ej)=zi ej−zj ei。
係数のz-degree dではd+2で割る。Cartan恒等式により全形式次数でδh=vであり、
全448個の別の固定多項式sign fixturesも独立照合した。一般性の根拠はCartan恒等式である。
この操作はS-linearでない。全780列の元の形式函数に個別適用してから、その列をS-linearに
並べる順序を全6240 vector nodeで確認した。normal form代表やscalar係数を先に外へ出さない。

## 完全128×972表示

Phi(s,b)=(I−nm)s+n d2b は ker Vgamma へ全射である。実際cycle vではmvがd1-cycleなので
mv=d2bとなり、Phi(v,b)=v。Phi(s,b)=0へmを掛けるとd2b=0、従ってb=d3cであり、さらに
s=n(ms)。ゆえに核は正確に(nS24,0)+(0,d3S168)で生成される。

各actual X列は(Xj,hj)へliftし、Phi(Xj,hj)=Xj。この780列を引いた後の残りの関係は
上記の核で全て表せる。従って保存された

    [[X,n,0],[h,0,d3]] : S972 → S128

はT_nil=ker Vgamma/im Xの完全表示である。Phiを掛けたfree44出力は[X,0,0]であり、
actual relationsを商にして初めて零となる。optional104×948の存在と、そのliteral complement
が未構成であることは区別する。全5408 scalar node、45 matrix node、6240 vector nodeの
型、source参照、全行列wiringとhashを照合した。

## 特殊化と長さ

元complexは各項がR0-freeなので、ell倍による短完全列のhomology長完全列から

    0 → N/ell N → T_nil → C'[ell] → 0

を得る。H2の未知性はこの箇所へ影響しない。受理済みD=QQ[[ell]]-finitenessより
dim N/ell N=mu_D N=L+3。C'のD-rank2とell商長4からdim C'[ell]=2。従って
length T_nil=L+5∈[16,22]。これは完全operatorの構成と長さの式であり、scalar成分、
minimal生成数、exact L、全nil生成元、nil三乗を計算したという主張ではない。

再現: `python code/ambient_B18_honest_nil_T_operator_review_20261104_v63.py`。
