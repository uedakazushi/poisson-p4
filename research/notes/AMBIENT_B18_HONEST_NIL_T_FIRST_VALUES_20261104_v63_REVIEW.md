# honest T_nil の最小生成数五：独立cache算術監査

2026-11-04。受理済み完全128×972表示の定数関係を、既存source cacheだけで監査する。
producerやその算術helperをimportせず、別の疎QQ monomial quotient算術を使う。
source solve、jet、点、GB、rank search、旧r/g quotientの再微分は行わない。

## 実source係数

C36=(k[z7,r]/m²)[tauA,tauB]/(tauA²−nuA,tauB²−nuB) と、B11上のC44を、保存された
monomial集合から一般的に再構成する。nu一次やC44のr²を落とさない。保存されたproper
kappa9、G8の最初五つ、実GL4を用いてdR9とD8×9を評価し、全44 Vgammaの一次係数と
第三sectorのr²係数を再現する。これらはcanonical Hadamard／r²除算の定数値そのもので、
元の正則quotientを微分する必要はない。m0の固定24minorの両側逆と行列式を検証する。

保存A12×9を新しくsolveせず、実N12について N A=[D;0] を全C36成分で確認する。
Nのorigin単元行列式と全12式 F12=A(dR9) も確認する。これによりAの値と向きは実source
から一意に定まる。Res4(A^t)をold B/Pの全640列へ適用し、P側の符号は変えない。
さらに全35 radical native syzygyを四C-basisで制限した140列を追加する。旧128列で
代用しない。全780列の定数・八方向一次係数を保存値に照合する。

## 定数消去と完全rank証明

固定18minorとその両側逆によりrank(X0)≥18、全780列の同じ18列による厳密factorization
によりrank(X0)≤18を確認する。各残余列は定数係数で選択18列を引いたものなので、
X0の像が厳密に零となり、(mX)_1=m0 X1 から未知m1 X0が相殺する。

元列ごとのHadamard／Euler contractionはS-linearではないが、QQ-linearである。
従ってこの定数列結合におけるh0は、m0 X1から読むconstant Koszul liftと一致する。
新Koszulの全84 wedgeの定数係数を抽出し、762残余それぞれの全24×8=192係数で
d2_linear h0=m0 X1を検証する。第三sectorでは最後のgeneratorがr²であるため、
最後成分のr一次係数が零であることも明示確認する。

W84の固定81minorについても両側逆と行列式を確認し、全762列のfactorizationで
rank≤81を証明する。source列を新しく選択していない。

## 結論の範囲

完全128×972表示のsection24単位対を消すと104targetsとなる。X0の18対を消した後、
W0の81対を消す。d3はminimalで定数関係を追加しない。従って

    mu_S(T_nil)=104−18−81=5.

これはS=R0/(ell)上、同値にR0上の最小生成数である。123単位対の消去後に5×849表示を
得られるが、その全形式scalar成分とliteral complementは構成していない。元128×972
表示は完全なまま保持する。T_nilの長さ、元nilidealの全生成元、nil三乗を結論しない。
rAとDVR parameter ellも区別する。

再現: `python code/ambient_B18_honest_nil_T_first_values_review_20261104_v63.py`。
