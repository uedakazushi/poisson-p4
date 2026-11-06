# 超平面＋既約四次曲面：停止時点の有限cofactor入口

2026-11-07。課題 `E_reduced_nonnormal`。
ユーザーの保存・停止指示に合わせ、本稿以降の新規計算は開始していない。
以下は既存源の確認と手導出した有限化の入口であり、全H+既約quartic支持の被覆ではない。
共鳴表の有限列挙、各型の全線形源、残Jacobiの全支持は未実施。
rootは有限化の方針を受理したが、以下の係数恒等式の独立算術証明書はまだ作成していない。

## 1. 全支持の正しい入口と既存四源

実退化式D=hF、hは単純な超平面因子、Fはhと異なる既約四次式とする。
h-Casimir liftを選び、h=1上で

    P=C+L(x)+Q(x),  Z=div P=c+Ax

と書く。C,L,Qは順に定数・一次・二次部分。
全PがPoissonであり、Zはそのmodular場なので[Z,P]=0、div Z=0。
従ってtr A=0。最高次QもPoissonである。

重要な点は、単純Hから `Q∧Z=0` は従わないことである。
重複Hの専用計算をこの層へ転用してはならない。

狭く確認した既知源は次の四つで、全支持のリストとは扱わない。

| 源 | 四次式Fの代表 | 性質・入口 |
|---|---|---|
| H_QUARTIC_ESSENTIAL5 | h⁴+h²x1x3+x1³x2−x1³x3+x1²x3²/2−x1²x3x4 | 本質五変数、全smooth23軌道。legacy/v06/nilpotent22.tex §quarticcomponent |
| MINIMAL_M | h⁴+3h²x²+12x³z | 本質三変数、全23成分。legacy/v07/minimal_thickness.tex |
| A0 | h⁴+h²xz−z⁴/4 | 本質三変数、全23成分。legacy/v07/nil22_v7.tex |
| DPM | h⁴+h²xz+p4(x,z) | 本質三変数、三つのcusp平面。legacy/v07/dplus_fullp4.tex |

いずれも既存証明は特定の全源又は軌道閉包の成分性を与える。
これらの存在から、任意のH+Fをβ=0 chartや同じ固定modular源へ移せるとは推論しない。

## 2. 無限遠のquartic Kは被約既約ではあり得ない

K=Pf(Q)=F|h=0と置く。Hが単純因子なのでK≠0。
QはC4上で一般に非退化な斉次二次Poissonテンソルである。
その逆形式ωは閉、weight零で、

    α=i_E ω

はP3上の非零の閉有理一形式となる。
分母はK、分子は三次式であり、Euler水平である。
α=0なら非退化ωに対してE=0となるため不可能。

Kが被約ならαはK上で高々単純極を持つ。
閉有理一形式の各既約極因子のresidueは定数であり、対数微分を引いた後に極は残らない。
P3には正則一形式がないので

    α=Σ_i λ_i dlog f_i,  K=Π_i f_i,
    Σ_i λ_i deg f_i=0.                                (1)

従ってKが被約既約quarticなら4λ=0、α=0となり矛盾する。
被約Kの可能な次数分割は1+3、2+2、1+1+2、1+1+1+1に限られる。

非被約K=Π f_i^(m_i)では、閉有理一形式の部分分数還元により、
(1)の対数部分にdRを加える。Rは同次次数零、極次数は各f_iで高々m_i−1。
したがってRの分母の総次数は4−Σdeg f_i≤3である。
この記述だけではQや低次lift C,Lの全支持は解決しない。

## 3. 非零半単純modular部の共鳴は有限本の直線へ絞れる

Affine4と全bivectorの非零scalar倍を用い、Zの半単純部を線形対角場

    S=Σ_i λ_i x_i∂i,  Σ_i λ_i=0

とする。残る冪零部はSと可換である。
[Z,P]=0から、有限次元の係数表現上のJordan分解により[S,P]=0も従う。
よってP_ijの単項式x_kx_l（x0=h、λ0=0）には

    λ_i+λ_j=λ_k+λ_l,  1≤i<j≤4, 0≤k≤l≤4             (2)

が必要である。
自明な恒等式を除き、trace零三次元空間上の有限なhyperplane arrangementを得る。
同じhyperplaneやtraceから冗長なものは一つと数える。
固有値一致λ_i=λ_jもこのarrangementに含まれる。

定数項又は一次項が非零なら、少なくとも

    λ_i+λ_j=0  又は λ_i+λ_j=λ_k                    (3)

という低次共鳴がある。
低次項が全て零ならP=Q、F=Kが既約quarticとなり§2に矛盾する。
従って(3)は必要である。

(3)の単一hyperplaneの一般点は、S4を除いて次の三型に尽くされる。
追加の(2)が成立する点はここでの「一般点」から除き、次節の交線に保持する。

### (a) pair-zero型 (a,−a,b,−b)

x,y,z,wの順としu=xy,v=zwと置く。全許容係数は

    Pxy=C+Au+Bv, Pzw=D+Eu+Fv,
    Pxz=p xz, Pxw=q xw, Pyz=r yz, Pyw=s yw.

交換項B,Eも保持する。
modular固有値のpair-zero条件はp+q+r+s=0。
v1=p+q、v2=p+rと置くと全Jacobiは次の五条件で尽くされる。

    BE=0, v2 C+BD=0, EC−v1 D=0,
    B(v2+F)=0, E(A−v1)=0.                            (4)

実Pfaffianは

    (C+Au+Bv)(D+Eu+Fv)+v1v2uv.                      (5)

B=E=0では、C又はDが非零なら(4)からv1v2=0で(5)は二次式の積。
C=D=0ならuvで割れる。
B≠0ならE=0,D=−v2 C/B,F=−v2。
C=0なら(5)はv因子を持ち、C≠0ならv1v2=0で再び二次式の積又は零。
E≠0の場合も対称にu因子又は二次式の積となる。
よって全境界を含め既約F支持は空。

### (b) zero-weight型 (0,a,b,−a−b)

零weightの座標をx、残りをy,z,wとする。
全許容係数はPxy=y(a1x+b1h)、Pxz=z(a2x+b2h)、Pxw=w(a3x+b3h)と、
Pyz,P yw,Pzwの三つのlogcanonical項だけである。
従ってPfaffianはyzwで割れ、Fは既約でない。

この一般点のzero eigenspaceは一次元なので、affine冪零部は零又はκ∂x。
後者を残しても[N,P]=0がx非依存を強制するだけであり、全Pfの因数分解は変わらない。
非零weight方向のtranslationはaffine座標変更で消せる。

### (c) triple-sum型 (a,b,a+b,−2a−2b)

全許容係数は六つのlogcanonical項と、P12の一次項κh x3だけである。
従ってPfaffianはx3x4で割れ、Fは既約でない。
この一般点では全固有値が非零で相異なるので、affine translationは消え、冪零部も零である。

以上より、非零Sを持つ既約F支持は、少なくとも一つの低次wallと、
それとは独立なもう一つの(2)のwallとの交線上に限られる。
trace零空間は三次元なので、二独立wallの交線をscalarで割った固有値比は有限個である。

## 4. 未実施の有限表と純冪零境界

次の計算入口は有限であり、未実施のまま保存する。

1. (2)をtrace零空間に射影し、primitive整数係数で同じhyperplaneを同一視する。
2. 一つの低次wallと全wallの交線を列挙し、零・同一wall・冗長なものを除く。
3. 得たλ四つ組をS4及び全scalarで同一視する。
4. 各重複固有値空間で全冪零Jordan分割を保持する。
5. zero eigenspace上でaffine translationをim Nを法として残す。
   Nのcentralizerはcoker N上で、零又は非零成分がある最大Jordan block長ごとの有限軌道を持つ。
   異なるblock長のtranslationを勝手に同一視しない。
6. 各固定Zについて全90係数のPにdiv P=Z、[Z,P]=0を課し、
   全Jacobiと実Pfaffianの既約性を最後まで保持する。

S=0でA≠0の場合はAの四つの非零nilpotent分割
(4)、(3,1)、(2,2)、(2,1,1)を扱う。
各分割でtranslation cの零／distinct block長のtopを分けると、順に2,3,2,3型、計十型。
これは源の型の数であり、成分数ではない。

A=0は除外できる。最高次Qのdivergenceが零であり、Pf(Q)≠0なら
Q#d Pf(Q)=0からd Pf(Q)=0となるが、非零四次斉次式には不可能。
従ってPf(Q)=0となり、Hが単純因子という条件に反する。
これはconstant modular場及び零modular場の両方を排除する。

## 5. 停止時点の射程

既存四成分の確認と§2–3の手による有限化が今回の到達点。
共鳴の全hyperplane数・有限λ代表数はまだ計算していない。
各係数型の全線形源・残Jacobi支持・既存成分との帰属も未実施である。
独立算術checkを行う場合は、まず(a)の全10係数Jacobi/Pfと、(b)(c)の全weight許容係数を
一つの小さい有限jobで確認する。次に初めて有限cofactor表を凍結する。
新しい源を既存のβ=0源だけに限定せず、全C,L,Qを保持すること。
