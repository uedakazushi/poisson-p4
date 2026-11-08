# Concurrent位数三の全Ore higher：strict四変数Poissonとbinary quarticへ

2026-11-08、PRIMARY G、producer純証明。独立review前。
新matrix、minor、GB、jet又はCASは無い。
whole sourceの正確な還元であり、G≤24又はE被覆の証明ではない。

## 1. 入口はfirst symbolだけがOre、higherは全て残す

base coordinates(a,b,c)、normal(u,v)、f=a b(a−b)とする。
rank二planeのnormal位数三、actualDは被約既約quintic、G3非零。
先行全first-coverのOre枝、すなわちraw first symbolが

`P0=∂c wedge W0(a,b,u,v)`

となるものを対象にする。
raw W0のa/b成分は(f_b,−f_a)、normal成分は
`(W0^u,W0^v)=M(a,b)*(u,v)^t`、Mは全linear matrix。
fullP1三十/fullP2九を保ち、初めからhigherをc-independentと仮定しない。

以下でy=(a,b,u,v)、E_yはそのEuler field、E=E_y+c∂cと書く。
δはδPのj成分がΣi∂iP^{ij}となるdivergence。
canonical `Pbar=P−E wedge δP/5` は全projective Poissonをstrict divergence-zero Poissonへ移す。
normalweightを保つこととBV identityによるstrict性は既canonical補題の通り。

## 2. canonical first symbolのnormal部分は非零

canonical P0は依然として `∂c wedge W(y)` であり、

`W=W0−(div_y W0/5)E_y`, div_y W=0。

W0のbase divergenceは零、div_y W0=tr Mである。
Wのnormal部分が零なら `M=(trM/5)I2`、traceを取りtrM=2trM/5なのでM=0。
従ってこの場合raw/canonical共にP0=Pf。

Pfのordinary conormal connectionは `β=(dlogf/3)I2`。
twisted Sym³ connectionのparallel条件は、各G3 coefficient gについて

`dg=(2/3)g dlogf`。

fの各prime line上でgのvaluationは2/3となり、nonzero polynomial gは存在しない。
normal位数三G3非零に反する。
従って本対象のcanonical Wは**少なくとも一つnormal成分が非零**。
このscalar first-symbol zero-leading補集合を隠れてinvertしていない。

## 3. 全canonical P1のc項とstrict weight一Jacobi

全P1はnormal degree一のbase-base項及びnormal degree二のmixed columnsである。
従ってcanonical P1を正確に

`P1bar=Q1(y)+c∂c wedge L(y)+c R(u,v)`

と書ける。
Rは `l1(u,v)∂a wedge ∂b`、l1はbinary linear。
Lはnormal-linear coefficientsを持つa/b方向vector field。
δ_y(cR)=0なのでdivergence-zero条件は

`L=−δ_y Q1`。

strict weight一Jacobiは `[∂c wedge W,P1bar]=0`。
その∂cを含まない三vector部分はsignを除き `W wedge R` だけである。
Wのnormal成分が非零、polynomial ringがdomainなのでR=0。
これはc²係数の推測やgeneric profile rankではなく、元全higherのliteral Schouten係数である。
ゆえに

`P1bar=Q1(y)−c∂c wedge δ_yQ1`。                       (1)

P2はnormal degree二のbase-base項だけで、初めからc非依存、δP2=0。

## 4. 明示Euler gaugeで全tensorをc非依存へ移す

V=δ_yQ1はnormal-linear coefficientsを持つa/b方向linear vector field。
(1)へactual Euler redundancy `E wedge V` を加えると、c依存項が正確に消え

`Q1(y)+E_y wedge δ_yQ1`

となる。P0/P2のnormalweightを変えない。
従って**任意の元full tensor**は同じprojective classのquadratic lift

`P=Q(y)+∂c wedge Wtilde(y)`                           (2)

を持つ。Wtildeも全quadratic vector fieldで、元sourceの全higher係数を保持する。
unique normalization又はsource上のfree actionを主張するものではない。
逆は(2)をそのままhomogeneous liftとして読むことであり、meromorphic座標変換やflatness仮定は無い。

projective JacobiをE wedge [P,P]=0で書く。
[P,P]のc係数は無く、E=E_y+c∂cなので、c係数から `[Q,Q]=0`。
残りは `E_y wedge [Wtilde,Q]=0`。
逆にこの二式なら(2)はprojective Poissonである。
従ってsource方程式は正確に

`[Q,Q]=0`, `E_y wedge L_Wtilde Q=0`。                  (3)

全350の一部分だけをPoissonと同値と呼んだものではなく、cの全係数比較による両包含。

## 5. strict Qはbinary二vector fieldsだけで記述出来る

元の全P1/P2とEuler補正のdegreeから、Qは必ず

```
Q=(a A1(u,v)+b B1(u,v)+C2(u,v))∂a wedge ∂b
  +∂a wedge U2(u,v)+∂b wedge V2(u,v),
```

A1,B1は全binary linear、C2は全binary quadratic、
U2,V2はnormal planeの全homogeneous quadratic vector fieldsである。
特にQ^{u,v}=0、mixed coefficientsはa,bを含まない。
Euler補正E_y wedge δ_yQ1もこのsource内に入る。

strict Jacobiの(a,b,u)/(a,b,v)係数を直接取ると

`[U2,V2]=A1 U2+B1 V2`。                               (4)

残二Jacobi triplesは自動零。
C2はこのstrict Q Jacobiへ寄与しないが、(3)のWtilde compatibilityとactualDへ残す。
全first symbolのf、normal connection、全leading targetはWtildeのweight零とactualD3へ戻す。
任意の(4)のsolutionが今回の元m三源を与えるとは言わず、(2)(3)と元のgraded条件が正確な全源である。

## 6. actualDはc一次、quartic coefficientはnormal binary

四変数orientation(a,b,u,v)でK4=Pf(Q)=−det(U2,V2)と置く。
(2)のactualDはoverall signを固定すれば

`D=c K4(u,v)+F5(a,b,u,v)`,
`F5=−[E_y wedge Q wedge Wtilde]_(∂a∧∂b∧∂u∧∂v)`。       (5)

係数二などを暗黙に落とさず、E wedge P²/2からの式である。
normal位数三ではG3はF5のnormal degree三部分で、c-independentとなる。
K4≡0ならDはcを含まない。この全補集合も残す。
K4≠0ならcommon point(a=b=u=v0)のlocal tangent coneはbinary quartic K4。
それだけからDがreduced又はirreducibleと推測しない。

## 7. K4非零の全root partitionは211/31/22/4

K=det(U2,V2)≠0とする。
(4)とLie derivative of U2 wedge V2から

`U2(K)=(div U2+B1)K`,
`V2(K)=(div V2−A1)K`。

従ってU2,V2はKの各distinct prime lineに接する。
homogeneous quadratic field T=(P2,Q2)のprojective componentはbinary cubic

`u Q2−v P2`。

Kに四distinct rootsがあれば、このcubicは四pointsで零なのでidentically零。
よってU2,V2はいずれもlinear multiple of E_normalとなり、det零に反する。
Kのdistinct rootsは高々三本。

三本の場合reduced cubic gがそのsupportであり、全quadratic logarithmic fieldsは

`T=l1(u,v)E_normal+t0 Jac(g)`。

projective cubicがgのconstant倍であることから直接従う。
二fieldのdetはg times linear form。
第四rootを持てないのでそのlinear formはgのいずれかのfactor。
従ってpartitionは正確に211。
二本なら31又は22、一本なら4。
これで非零quarticの全補集合を尽くし、four-distinct-rootを勝手にgenericと仮定しない。

本稿はこれらrepeated binary strataのwhole Ore/Wtilde compatibilityをまだ分類していない。
既m四normal-leading分類と同一のsourceであるとは主張しない。

## 8. Free-triple flatness符号の独立receipt

関連sourceの符号を独立に記録する。β=[[p,0],[r,q]]のrow作用ではcurvatureは

`dβ−β wedge β`, lower entry `dr−(q−p) wedge r`。

本producerの先行root geometryが用いた `dr+(p−q) wedge r` と同じ式である。
tracefree Zに対しω=i_E i_Zvol、r=−ω/(3f)とすると

```
dω=3i_Zvol,
dr=−i_Zvol/f−dlogf wedge r.
```

θ=p−q−dlogfはEuler値−3を持ち、

`θ wedge r=i_Zvol/f+θ(Z)i_Evol/(3f)`。

従ってflatnessはθ(Z)=0。
triple q=2p−dlog(f/h)ならθ=−p−dlogh、すなわち

`(p+dlogh)(Z)=0`。

先行concurrent Tj/central純証明のχ条件と一致する。
今回のcorrected符号から導いたので、別sourceの旧誤ったplus-sign分類を依存にしない。
full non-Ore算術も全元P0 Jacobiを直接検算している。

## 9. 正確な次の義務

全concurrent m三残余は、未列挙c-independent SF有限first coversと
全repeated Ore first symbolsを、(2)(3)(4)(5)のactual sourceへ接続したもの。
K4零と全211/31/22/4を保ち、Wtildeの元f/leading targets/全higherを戻す必要がある。
本稿だけでは全fixed-source次元、PGL像次元、componenthood又はE所属は得ていない。
同じprofile/minor/GB/jetを追加するのではなく、このstrict binary sourceとactualDの構造を次手とする。
