# generic B18のactual odd ideals：discriminantで割った厳密な形式型

2026-11-03。CA/CB/Nのsource critical schemeの独立受理済み完成型を使い、
actual imageのodd coefficient idealを**各discriminantを法として**分類する。
二つのCA/CB座標変換を同じtarget上で整列させた結果ではなく、actual Dの完成型は未解決である。

## 1. Nには別の28次元cone parameter baseを作る

CA/CBに共通な底は29次元cone底R29であり、source次元28のNに流用できない。
そこでN自身の旗9＋fiber19というsmooth28次元sourceの微分を全150係数で計算する。
source順序をfiber19列、次いで既存の旗9列とする。
列22（旗E_xh）のみを除いた27列は、次の実target27行でdet−1152/5を持つ。

    5,15,16,20,22,23,31,35,50,65,105,106,108,109,
    110,112,113,117,118,119,120,121,123,125,127,132,140.

M=E_xh−2E_hzとし、LP=Lie_M(P0)とする。
LPは固定N-fiber19空間に入り、旗Mとfiber速度−LPからなるkernelは列22成分1である。
実有限source曲線

    exp(pM)_*(P0−p LP)

の一次target項は零、二次項はW=−Lie_M²(P0)/2である。
上の27行を固定する陰関数補正後、実targetの行27の二次係数は**−3**となる。
したがって27行と行27の計28実target座標でsourceの閉fiberを切ると、
完成局所環は厳密にC[[p]]/(p²)となる。
高次項は非零二次係数を持つ単元倍なので長さを増やさない。

完備有限性と正則sourceのCM性から、N-normalizerはこの別のregular R_N,28上有限自由rank2である。
行106は選択27行に含まれ、中心で値1なのでprojective版では一つを除きR_N,27となる。
このN底の構成はCA/CB共通底との同一視を含まない。

## 2. quadratic normalizerのdifferentとactual odd ideal

三族を一つずつ考える。CA/CBには既証の共通R29、Nには上記のR_N,28を使う。
normalizer BはR自由rank2で、閉fiber basisを持ち上げた1,pを基底に取れる。
traceの半分を引いてt=p−Tr(p)/2とすれば

    B=R[t]/(t²−nu),  actual image A=R⊕I*t

となる。Iは残るactual target出力の全odd係数が生成するRのidealである。
底Rはactual target座標から作られているのでR⊂Aであり、任意の抽象normalizer底をimageへ混ぜていない。

Bは正則でtが中心で一次のsource parameterである。
上のhypersurface表示でBの正則性を見ると、nuのRでの一次項は非零でなければならない。
したがってnuはRの正則parameter、tもBの正則parameterである。

Rに相対的な微分はOmega_(B/R)=B/(2t)*dtで、Aの追加元i*tの微分はi*dtである。
従って全normalization射のramification Fitting idealは厳密に

    Icrit=(t,I)⊂B,
    B/Icrit = R/(nu,I).

tは非零一次項を持つので、t∉m_B*Icritである。
ゆえにtはIcritのminimal generatorの一つとして選べる。

## 3. modulo discriminantでの全formal ideal型

ここでIbar=(I+(nu))/(nu)⊂Rbar=R/(nu)とおく。
これは**idealの像**であり、module base change I/nu*Iと同一視しない。
I∩nuR=nuIを仮定しておらず、Aのnu=0 fiberの全ring型も推論しない。

独立受理済みのsource critical完成環は、coneで次のとおりである。

    CA,CB: C[[g19,w,u,v]]/(wu,wv,uv,u³),
    N:     C[[g19,w,v]]/(wv,v²).

一方、RbarはCA/CBでregular28次元、Nでregular27次元である。
上のcritical商の埋め込み次元はそれぞれ22と21なので、六つの独立一次消去を行える。
最小Cohen presentationの一意性（cotangent基底のliftと形式逆写像）により、
Rbar内の**全Ibar**は独立な形式ambient座標変換の後で次のidealとなる。

    CA,CB:
      Rbar ≅ C[[g19,s1,...,s6,w,u,v]],
      Ibar = (s1,...,s6,wu,wv,uv,u³).

    N:
      Rbar ≅ C[[g19,s1,...,s6,w,v]],
      Ibar = (s1,...,s6,wv,v²).

これは単にcritical商の抽象同型を記しただけでなく、regular ambientの次元を含めた
ideal pair (Rbar,Ibar)の形式同値である。
ただしCA/CBの二つの同値変換は**独立**であり、共通R上の相対位置は決定していない。

## 4. actual odd idealに関する有限な帰結

上の全生成元はminimalであり、mu(Ibar)はCA/CBで10、Nで8となる。
同値に、source critical idealのminimal生成数11,11,9から、
上記のregular tを一つ除いた数である。
I→Ibarの全射により、full odd ideals自身のminimal生成数もそれぞれ10,10,8以上である。
したがってCA/CBで九個、Nで七個のregular odd generatorsを選ぶ旧Ac型のshortcutは成立しない。

この結果はdiscriminantを法とした全ideal型とminimal生成数の下界を与える。
nu方向のdeformation、二つのimage kernelの和、actual Dのnilradical・associated primes・depthは依然別の問題である。

* code: `code/ambient_B18_odd_discriminant_reduction_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63/`

Nの全150×28微分、27単元minor、kernel、実有限curve、全150二次補正を保存する。
10/10/8のminimal生成数は別のMacaulay2計算でも確認する。
