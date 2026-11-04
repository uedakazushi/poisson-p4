# 固定 F32003 多項式による R434 と Saux の除法

2026-11-04。新しいparameter点やprimeを使わない。一つの多項式環F32003[a,e,d]で、
受理済みR434と7022項Sauxの全QQ係数を還元する。全QQ分母のp単位性、Rbarの非零と
全次数20保存、Sauxbarの全次数61保存を最初に検査する。因数分解はこのpacketで行わない。

graded reverse lexicographic orderで一回だけ `Sauxbar=Q*Rbar+rem` を計算する。
Qの全次数は41以下、remの全次数は61以下で、remのどの単項式もRbarのleading monomial
で割れないことを検証する。積は全係数で掛け戻す。rem≠0ならRbarはSauxbarを割らない。

QQからの帰結は点試験ではなくGaussの補題である。Rにはp単位の係数があるため
Z_(p)[a,e,d]でprimitiveであり、QQ上のR|Sauxはp整の商を持つ。従ってmod pの
非整除はQQの非整除を証明する。さらに別の独立R既約性gateが受理された場合にだけ
gcd_QQ(R,Saux)=1を結論する。多変数ideal(R,Saux)=1や空の交叉は結論しない。

支持予測はQ≤C(44,3)=13244、全次数61dense≤C(64,3)=41664。remは一つの
degree20 leading monomialの倍数を含まないので≤41664−13244=28420。
remをstream保存して、T=Sauxbar−remを作ってからremの係数memoryを解放する。
T≤35442、掛戻しRbar*Q≤41664。input7456、Q13244、T35442、product41664、
output-sized scratch41664、reserve10000を含む149470項が最大上界となる。
division内部支持はinstrumentしないので、その不確実性と60秒/2GB制限を併記する。
全係数は16bit、全資料32MiB以下。Python親もnative実行前にpoly配列を解放する。

一回の固定除法だけを実行する。非零remainderの際も、他のminor、Saux factor、
prime、QQ除法、一般gcdへ自動的に進まない。結果が零又はguard停止ならそのまま保存する。
この判定は実pure pivot chartの補集合についての代数的情報に限り、残る5×12のrank、
全FULLU57、global exhaustionを結論しない。
