# 三軸 conductor の完成局所環と二次函数の剛性

2026-11-03。独立した概念的監査。以下は**仮定を明示した抽象的な完成局所環の定理**である。数値計算・GBは行わず、CA/CBへのgroup product・cover embedding・conductor pushoutの適用は別証明とする。

## 定理

Cを剰余体が複素数体である1次元Noetherian完備局所C-algebra、極大idealをmとする。次を仮定する。

1. CはCohen–Macaulayで、embedding dimensionは3。
2. C_redは三つの滑らかな1次元枝を持ち、その三つの接直線は独立。
3. 三つのminimal primeにおけるCの局所長は、二つが1、一つが2。
4. graded algebraとして

\[
\operatorname{gr}_{\mathfrak m}C\simeq
\mathbb C[w,r,s]/(wr,ws,sr^2).
\]

このとき、完成局所C-algebraとして

\[
C\simeq M:=\mathbb C[[w,r,s]]/(wr,ws,sr^2).
\]

長2の枝はs軸、単純二枝はw,r軸に対応させられる。以下の構成では、最初のgraded同定に対し線形項が恒等な同型を選べる。

さらにq∈m²が二つの単純な被約枝上で恒等的に零であり、その**初期二次類**が

\[
q\bmod\mathfrak m^3=c(s^2-rs),\qquad c\ne0
\]

なら、枝の標識を保つ形式的座標変換により、環と函数のpairは

\[
(C,q)\simeq(M,s(s-r))
\]

となる。qの単位倍への置換だけを主張するのではなく、qそのものをこの函数にできる。

## 1. CM性から冪零部分の作用を決める

N=Nil(C)とし、長2のminimal primeをP_Sとする。単純二枝のgeneric点ではCは体なのでNは零、P_SではCの局所化は長2のArtin局所環で、その極大ideal=N_(P_S)は1次元かつ二乗零である。

1次元CM局所環Cは有限長の非零部分加群を持たない。実際、そのような部分加群があればmがCのassociated primeとなりdepthが0になる。従って、Cの部分加群で全minimal primeにおいて零となるものは零である。

任意のa∈P_Sについて、aNはP_Sで零になる。長2のArtin局所環ではaの像が極大idealに入り、極大idealの二乗が零だからである。他のminimal primeではN自体が零である。よってaNは有限支持を持つCの部分加群となり、aN=0。従って

\[
P_SN=0,\qquad N^2=0.
\]

このためNの作用はC/P_Sを経由する。S枝が滑らかなのでC/P_S≅C[[s]]はDVRである。Nはその上の有限生成moduleで、s-torsionがあればCの有限長部分加群になってしまう。ゆえにNはtorsion-free、generic rankは1であり、

\[
N\simeq\mathbb C[[s]]\varepsilon
\]

と自由になる。C[[s]]のsectionを先に仮定する必要はない。S枝のparameterの任意のliftを取り、完備性からC[[s]]→Cを定めれば、C→C/P_Sとの合成が恒等なのでsectionが得られる。

## 2. 二重の枝の軸と reduced 三枝の座標

graded ringの被約支持は三つの座標軸である。滑らかな被約三枝の接直線はこの支持に含まれ、独立なのでそれぞれが三軸に対応する。

長2の枝がs軸に対応することも、単に標識として仮定しなくてよい。Nの非零元εの初期類は冪零である。上のgraded ringのnilradicalは `(rs)` であり、次数dの非零な冪零類は `c rs^(d−1)` である。その線形annihilatorは正確にspan(w,r)である。一方P_SN=0なので、S枝に垂直な二つの線形方程式はこの初期類を消す。従ってSの接直線はs軸である。残る単純二枝をw,r軸に沿って命名する。

embedding dimension3の最小正則表示 `C[[w,r,s]]→C` を取り、その被約三枝のparameter表示をγ_W(w),γ_R(r),γ_S(s)とする。接ベクトルをそれぞれ標準基底に正規化すれば

\[
\Phi(w,r,s)=\gamma_W(w)+\gamma_R(r)+\gamma_S(s)
\]

のJacobianは原点で恒等である。形式逆函数定理によりΦはambientの自己同型であり、三つの座標軸を実際の三枝へ送る。従って線形項を変えずに

\[
C_{red}\simeq\mathbb C[[w,r,s]]/(wr,ws,rs)
\]

とできる。この変更はassociated gradedの指定した表示を保つ。特にNはC内の `(wr,ws,rs)` として生成され、N⊂m²である。

## 3. rs がnilmoduleの生成元である理由

この座標でP_Sの像はw,r,Nを含むので、wN=rN=0である。DVR上のNの生成元εを一つ取ると

\[
rs=\gamma(s)\varepsilon.
\]

もしγ(0)=0なら `rs∈sN⊂m³` となる。しかし指定されたgraded ringではrsの次数2類は非零である。従ってγは単位であり、**ε=rsを生成元に選べる**。

ここで `in(N)=Nil(gr C)` やnilradicalのbase changeを仮定していない。NのDVR自由性、N⊂m²、rsの非零初期類だけを使っている。

wr,wsはいずれもNに入り、`wr=α(s)ε`、`ws=β(s)ε` と書ける。結合律とwε=0により

\[
s\,wr=w(rs)=w\varepsilon=0.
\]

sはN上単射だからwr=0、α=0となる。さらにgrの次数2ではws=0、rs≠0なのでβ(0)=0である。従ってβ=s b(s)と書ける。

\[
w'=w-\beta(s)r
\]

と変更すると、線形項は変わらない。そして

\[
w's=ws-\beta rs=0,\qquad
w'r=wr-\beta r^2=-b(s)r^2s=0,
\]

最後の零はrε=r²s=0から従う。従ってw',r,sはCを位相的に生成し、monomial ring MからCへの全射を与える。この全射はassociated graded上で恒等同型なので、核は零である。非零な核元があればその最初の非零初期類が消え、graded同型に矛盾するからである。これで環の剛性を証明した。

## 4. 函数qの完全な正規化

Mの全元は一意に

\[
c+w f(w)+r g(r)+s h(s)+rs j(s)
\]

と書ける。qが二つの単純枝で零ならc=f=g=0である。q∈m²を使うと

\[
q=s^2u(s)+rs v(s).
\]

初期二次類からu(0)=c、v(0)=−cで、両方とも単位である。複素数体上でuの形式的平方根A(s)を取り、

\[
S=sA(s),\qquad R=-r\,v(s)/A(s),\qquad W=w
\]

とする。A(0)≠0、v(0)≠0なのでambientのJacobianは可逆である。`WR=WS=SR²=0` は元の三関係から従い、逆変換も同じ形なのでMの自己同型である。そして厳密に

\[
S(S-R)=s^2u(s)+rs v(s)=q.
\]

この変換は三つの被約軸をそれぞれ保つ。qの定数倍も残らない。

付随する整合性として、U=S、V=R−Sとすればq=−UVであり、

\[
C/(q)\simeq\mathbb C[[W,U,V]]/(WU,WV,UV,U^3).
\]

これは抽象pairの帰結であり、実際のCA/CB coverと同定したという意味ではない。

## 5. 仮定を弱める場合の注意と反例

* q₂は **q∈m²の下での初期類** を意味しなければならない。単に二次Taylor部分だけを指定するなら `q=s+c(s²−rs)` は二つの単純枝で零だが、order1であるためorder2の `s(s−r)` には自己同型で移らない。
* associated gradedの表示だけでは三つの被約枝とgeneric lengthsは決まらない。例えば

\[
\mathbb C[[w,r,s]]/(wr,ws,sr(r-s^2))
\]

は同じassociated gradedを持つが、w軸と、平面内のs=0、r=0、r=s²という四つの被約generic branchを持つ。平面の主idealの初期式はsr²であり、wはr,sを消すので、上のgraded同一性が直接分かる。これは定理の三枝・generic lengths仮定を満たさず、その必要性を示す。
* 一般の高次元底上のrank1 torsion-free moduleはDVR上のように自由とは限らない。この証明をg19を残したまま無条件で適用しない。実際のconductorとcoverをgroup作用により形式的にproductへ分解する段階は別に必要である。

## 範囲と解析的意味

定理は完備局所C-algebraの形式同型を与える。出発点が収束複素解析germで、同じ滑らかな三枝とgeneric lengthsの仮定を満たす場合は、枝の直線化、DVR module、単位の平方根を収束環C{s}で同様に構成でき、収束する解析同型も得られる。任意の形式座標変換の収束を無根拠に仮定してはいない。

この証明書は、CA/CBの実conductorがこれらの仮定を満たすこと、g19とのproduct、tracefree coverの標識、normalization内のembedding、conductor pushoutの正確なambient像を一括して受理するものではない。それらの適用を確認した後に初めて、実際のcompleted imageの型へ用いる。
