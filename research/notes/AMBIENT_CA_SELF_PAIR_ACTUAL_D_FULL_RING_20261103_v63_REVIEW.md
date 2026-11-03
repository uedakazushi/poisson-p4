# CA自己旗族に沿うactual D全環の独立承認

2026-11-03。対象は `AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63`。
全proper fiber・全81 fiber-product kernelを含むmeeting-double独立監査と、
QQ(a)の全150真二次差を含むD-twojet独立監査を入力として固定する。
その物理計算は再実行せず、両ideal包含、同時座標、全ambient上界、nil/Ass/depthと大域化を監査する。

## 同時座標の正当性

CAの二smooth29枝をB1,B2、meetingをdouble smooth20のCとする。
独立入力は `B1 intersect CB=N` がsmooth28、`B2 intersect CB=S` がsmooth20という
**scheme-theoretic**な等号を与える。従って `C intersect N=S` もschemeとして成立する。
第一枝内のNの方程式はC上でunit*rhoとなる。これを第一枝のrhoとし、同じCのclassを
第二枝へliftすれば、N=(rho)とmeeting kernel `(u8,rho²),(v8,rho²)` を同時に保てる。
単に三つの接空間の次元を合わせて座標を仮定していない。

本packetのtauは `(0,rho²)`。第一枝はv=tau=0、第二枝はu=0,tau=rho²である。
旧meeting監査のtau_old=(rho²,0)からは `tau_old=rho²-tau` となる。
独立の整数word-polynomial engineで、この置換が全81関係を符号以外同じ関係へ移すこと、
および全81関係が両枝で零となることを確認する。

## 両ideal包含と追加ambient式

CBはsmooth29、Nはそのsmooth divisor、SはN内のsmooth20なので、
`B=C[[s20,u8,f]]`、N=(f)、S=(f,u8)と取れる。
CAのrho,v8,tauのliftはNで零なので、CB上ではrho=fd、v=fb、tau=fcと書ける。
u8の像はN上のSの正則座標であり、fと共にBの正則座標をなす。

d tauは二CA枝の接空間の和37を殺すが、whole CA tangent38の追加方向を検出する。
受理済みの全S上の `T_CB subset T_CA` とrank triple38から、CB方向でd tauは非零。
従ってc|S≠0、cはlocal unitである。
minimal CAの38 generatorを元ambientへliftして使える。追加のambient共通式を捨てる必要はない。

全81関係中の `u_i*tau` と `tau*(tau-rho²)` のCB像は

    c*f*u_i,  f²*c*(c-f*d²)

であり、cとc-fd²はunit。従って `f*(f,u8) subset Q` を得る。
独立コードは全81像を自前の可換word算術で再構成し、保存値と全係数を照合する。

逆包含は最小81式だけの検査ではない。
N⊂Dから、**任意のactual CA ambient equation** のCB像はF=f*kと書ける。
稠密開Sの全点でT_D=T_CBだから、`dF|S=k|S*df=0`。
dfは非零、Sは被約なのでk|S=0、従ってk∈I_S。
全ての追加ambient equationも含めて `Q subset f*I_S` となる。
従ってQ=fI_Sは全idealの等号であり、同じ有限変数表示を仮定した循環はない。

## nil・associated prime・depthと次数

全完成環は `B/(f²,f*u8)`。
monomial LCMによって `(f²,f*u8)=(f) intersect (f²,u8)` を独立確認する。
第一成分はprime、第二成分はS-primary、両者はirredundantなのでAssは正確にNとembedded S。
nilradicalはfB/fI_S≅O_S、平方零、S上のnil moduleのgeneric長は1。
embedded primary quotient自体のArtin長と、このnil moduleの長を混同しない。

g20を除いたcoreの標準monomialは任意のu8 monomialと単独のfだけ。
Hilbert seriesは `1/(1-T)^28+T/(1-T)^20`、dim28、edim29、e1。
coreではfが全maximal idealに殺される非零元なのでdepth0、g20はregularなのでdepth20。
projectiveではscalarのsmooth factorを除いてdim27/edim28/depth19/e1。
e1をreduced/CM判定に使っていない。

## embedded成分の大域的閉包

Sはactual G×a像の代数的閉包のsmooth開部分であり、NとI_Sは代数的に定義される。
局所CB内のfもNの代数的Cartier方程式に取れる。
証明された完成ideal等号 `Q_hat=f*I_S_hat` は、忠実平坦なcompletionにより
元の代数的local ideal等号へdescendする。
従ってnilのannihilator I_Sは元のlocal actual Dにもassociated primeとして存在する。
Noetherian schemeのAssのlocalizationから、その閉包はglobal actual Dのassociated component。
SはN内のproper20次元subsetなのでembeddedであり、cone20/projective19である。
単に未代数化のformal associated primeを大域化したわけではない。

承認範囲は非零parameterのある非空開とそのG飽和、かつverified meeting-double unit open内。
全a≠0、a=0/B18原点のD、他のglobal primary component、CB自己旗族での別成分は主張しない。
旧R19は別のminimal componentであり、このembedded Sと同一視しない。

独立コード・全hash・証拠は
`code/ambient_CA_self_pair_actual_D_full_ring_review_20261103_v63.py` と
`reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/` に保存する。
