# 全PGL二不変量と二次の挟み込みの独立監査

2026-11-04。対象は凍結された `Q_NEW_FAMILY_GLOBAL_TWO_INVARIANT_PREFLIGHT_20261104_v63`。
対象の体は実際の25次元Poisson成分の **full PGL5 invariant field** であり、
15次元firstjet像のSp4商体とは区別する。producerを呼ばず、既存全affine tensorから
有限modular flow、H上の全二次係数、log形式、adjugate pairingsを再構築する。

## 評価点とscalarの除去

全tensorのdivergenceから直接

    dy/dt=4u/3,  dz/dt=u*x+u*v*y+1/48,  dx/dt=dw/dt=0

を得る。指数3の有限flowについて、全六成分で
`P(flow(x))=Dflow*P(x)*Dflow^t` を確認した。`t=-3y/(4u)` により任意の一般
affine評価点をy=0へ移せ、u,vを変えない。y=0上の六sourceではzとBが
非零係数−4u/3の一次変換で互いに回復される。全valueと任意scalar cでの
Darboux frameも直接照合する。従ってこのsourceは特定の評価点だけのものではない。

受理済みの実際のtarget逆写像・指数1からrはtarget上のSp4不変な有理関数である。
そのsource pullbackは受理済みsource derivationの定数である。定数体を全て生成する
有理section上では、受理済み圧縮の逆式により `r=u²/1600`。またδu=0なので、
この等式は全六sourceで成り立つ。cも自由だからtensor scalarを除去する。
任意のPGL変換後の固定点での計算は元tensorの別の評価点での計算となる。
valueとperspectiveの正常化間の残る差はSp4なので、`A=1600r=u²` は実際のfull
componentの有理PGL不変量である。全境界でregularとの主張ではない。

## intrinsicなfoliationと第二不変量

一般のquinticは3H+QでHは唯一の三重hyperplane、Qはsmooth quadric。
Q|Hの唯一のvertexからの射影はP2(x,y,w)をintrinsicに定める。
全二次係数を独立抽出して `P|H=∂z∧V`、Vのz非依存性を確認した。
従ってcharacteristic distributionの射影はintrinsicなline fieldである。

`l±=x±8w`、`q=-x²+64w²+64uxy+32uvy²` とすると、独立計算した
log分子はEuler方向とVをともに消し、非零である。log形式のresidueは
`1,u/16−1,−u/16−1`。他の幾何的generic invariant lineへの制限が零なら、
各点の整数位数A,Bは `A+(u/16)B=0` を満たす。uは超越的なのでA=B=0。
二つの比が定数となり、線はO=(l+=l−=0)を通る一方、q(O)=32uv≠0が
第二の比の定数性に反する。よって一般の幾何的函数体では二本だけである。
特殊共鳴parameterについて一律の不存在を述べない。

authorのinverseとは別のadjugate法でpairingsの分子を得た。

    self+=self−=−65536u²,  cross=4096u(v+16u).

従って `J=256 cross²/(self+ self−)=(v/u+16)²`。
この式はconicと各lineのscalar倍、二本の交換、PGL3座標変換に不変である。
幾何的genericで一意なunordered pairに対し対称だから、QQ上にも有理下降する。
これはfull tensorの葉層からの構成で、単なるVのsampleに対する値ではない。

## 体の結論と範囲

全tensorの既知の同時符号変換も全六成分で再照合した。dense full sourceの
pullbackとこの変換により、k=QQ,Cに対し

    k(A,J) ⊂ E_full ⊂ k(A,b),  A=u², b=uv,
    b²+32Ab+(256−J)A²=0.

判別式4A²JのJ-adic位数は1なので右端/左端は正確に二次。
中間体は両端点しかなく、どちらも純有理二変数体である。
従って**実際の一般full PGL moduli体はQQ/C上で純有理・超越次数2**。
X25自体の有理性はこの論証に使わない。

この凍結監査はどちらの端点か、追加sigmaの存在、全特殊parameterや有限isotropy群を
決定しない。後続のsigma証明は別packetとする。新しい点・prime・高次数探索・GBはない。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_global_two_invariant_review_20261104_v63.py`
（60秒/2GB上限）。
