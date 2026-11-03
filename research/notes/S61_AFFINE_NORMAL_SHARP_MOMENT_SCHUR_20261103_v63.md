# Sharp532 controlled localizationと135本のscalar algebra表示

2026-11-03。先の536構成は正しいまま保持し、不要な(h−a)を除いて改善する。新delta0では、**一つのgeneric a、一つの非零因子、670本degree≤10または135本degree≤34のscalar理想**による正確な未知残余の空性判定を得る。unit自体はまだ証明していない。

## 1. Monic kernelから自動的な四列の選択

前と同じr(a)=sum a^i r_iを使うが、四列を

    b_j=T_h^j*r(a), j=0,1,2,3

とする。fixed bottom四行のdetをdelta0(a)とする。未知点pでr(a)がideal(h−h0)を生成すると、そのmultiplication kernelはmonic degree4のg4=G/(h−h0)のspanである。従ってdegree≤3のdomainは**常に**そのkernelの補空間であり、(h−a)による別のcomplementは不要である。

r_i=(h−h0)s_i、D=k[h]/g4とおくと

    delta0(a)=Norm_D(sum_i a^i*s_i)

が正確に成立する。s_iはDのunit idealを生成するので、このnormはaについて非零、次数≤4*133=532である。Dの非被約性は妨げにならない。旧constructionにあった追加g4(a)の因子を除いたので、bad値上界は536から532へ下がる。

producerはsimple、多重根、G=h^5/r0=h²/r133=hの三つのexact fixtureで全5×5 multiplication matrixとnorm identityを検査する。最後はdelta0=a^532となり次数上界を達成する。

## 2. 非零因子と正確な一つのunit criterion

全未知幾何点でdelta0(a,p)が非零多項式なので、その係数idealはBのunit idealである。Ass/McCoyによりdelta0はB[a]とB(a)の非零因子で、B→B(a)_delta0は単射である。全associated componentsとnilpotentsが残り、delta0を反転してからunitを証明すれば**B=0が正確に従う**。

known-five factorでは全Phi entriesが零なのでdelta0=0であり、A(a)_delta0=B(a)_delta0となる。これは未知componentを丸ごと落とさないことが証明されたcontrolled localizationである。全閉点を一chartが覆うとは主張しない。

全点coverが必要ならdim B≤2と三つの独立generic aで三chartを得る。また0,...,532の533個の固定QQ constantsなら、dimension boundなしに全点を覆う。旧537固定値によるcoverも正しいが、新構成では533で足りる。

## 3. 670本degree≤10の表示

四b列をpivotとしたSchur式

    det[b0,b1,b2,b3,Phi_column_l], l=0..669

を用いる。delta0を反転した環でこれらのscalar idealは実A(a)_delta0を正確に表す。b列のformal parameter次数≤2なのでdelta0≤8、各式≤10である。a次数は≤532である。

## 4. 135本のalgebra表示

M_r=r(a)(T_h)の五列全てを使い、N=det M_r=Norm_C(r(a))とする。初めの四列のbottom block Pが可逆なので、C/(r(a))ではbasis1,h,h²,h³,h⁴の後ろ四元を1の倍数へ解ける。特に

    h=tau=U/delta0,
    U=−(toprow_b*adj(P))_0.

fifth relationはN/delta0である。したがって**unital algebraとして**

    C/(r(a)) localizeddelta0 = R'_delta0/(N)

となる。module cyclicityだけを自由性と混同していない。1が実際のunit generatorであり、全commutative monic algebraのquotientであることを使う。

r(a)は元r_iの結合なので、実Aはこの環に全134式r_i(tau)を課したものである。r_iのh次数≤4だから、delta0^4を掛けて

    N,
    sum_(k=0)^4 r_i,k*U^k*delta0^(4-k), i=0..133

という**135本のscalar generators**を得る。Nのformal次数≤10、Uとdelta0≤8、残る134本≤34である。式数を減らす代わりに次数上界が大きくなるtradeoffがあり、670本degree10の選択肢も残す。

producerはcoefficient-ring fixture G=h^5+c*h+b、r=h−xでdelta0=1、tau=x、N=−(x^5+c*x+b)を正確に確認する。これはalgebra表示のimplementation checkであり、元S61のunitを示す標本ではない。

## 保存と残件

`code/s61_affine_normal_sharp_moment_schur_20261103_v63.py` と `reports/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63/`。全134source・global quintic・既証明normal閉埋込みをpinする。巨大delta0や135個のnormal展開、actual scalar unit、全w飽和は未計算である。modp unitやnilpotenceだけから一般QQへ移す推論はしない。
