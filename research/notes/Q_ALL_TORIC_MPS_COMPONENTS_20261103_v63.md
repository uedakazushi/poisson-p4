# 全 Q_a を通る三十八の toric MPS 成分

2026-11-03。全有限複素数 a に対して、canonical 二次 Poisson tensor

    (Q_a)14=a*x0²,  (Q_a)23=−2*x0²/9,
    (Q_a)24=4*x0*x1/3,  (Q_a)34=2*x0*x2/3+x1²

の射影類は、既存 MPS カタログの **T01–T38 の全三十八成分**に属する。
CA、CB への既証の所属と合わせ、全 Q_a を少なくとも四十個の相異なる既存成分が通る。
a=0 の rank2 tensor も含む。係数 −2/9 が固定されているので Q_a は全 a で非零であり、
a↦[Q_a] は全 A1 上の射である。

これは正の所属定理である。通過成分・局所枝の本数を尽くした定理ではなく、
T39、T40 や他の成分への非所属を主張しない。全 Poisson 成分の既存下界57は増えない。
また Q_a の特殊一次ジェット上の結論であって、一般 beta の affine fiber が
既知五点だけであるという未解決問題を解いたものではない。

## 実際の canonical arc

rho を変数として

    s=3(rho−1)/2,  a=8rho/[9(rho²−1)],
    mu=−4/[9s²(s+3)]

と置く。regular nilpotent N_s の四つの subdiagonal は
1,s,(s+3)/2,−9/(s+3) である。互いに異なる d=(0,1,z,x,y) に対して
N_t=N_s+t*diag(d)、W_k=N_t^k−tr(N_t^k)I/5 とし、線形場の wedge を用いて

    Pi_t=mu*(alpha*t² W1∧W2 + beta*t W1∧W3
             +rho W1∧W4 + W2∧W3)

を取る。行列冪は可換かつ tracefree なので、これは actual canonical Poisson tensor の
t に関する多項式族である。追加二項は t=0 で消え、全150係数で Pi_0=Q_a となる。

t≠0 では単位下三角 eigenframe により対角化できる。
共通非零 scalar mu*t5 を除いた対数型 bivector の六係数は

    gamma_ij=alpha(di*dj²−di²*dj)+beta(di*dj³−di³*dj)
             +rho(di*dj⁴−di⁴*dj)+(di²*dj³−di³*dj²) .

period は gamma 自身でなく、その線形 Pfaffian dual

    (gamma34,−gamma24,gamma23,gamma14,−gamma13,gamma12)

である。Pf(gamma)=rho*prod(di)*prod(dj−di) なので admissible 開で rank4 である。
標識置換後の period が対象の保存 period 線形空間へ入れば、Pi_t はその actual toric-boundary
族の PGL 変換である。この族は対象成分に含まれ、period の別の開条件が落ちる点も
Pfaffian-dual 多項式族の閉包として同じ成分に入る。

## 全三十八成分を覆う証明の分割

| 対象 | 証明された実 source |
|---|---|
| T01 | 可換行列冪による基本 toric arc |
| T02,T03,T16,T17 | alpha,beta,rho と固有値の明示的な有理函数族 |
| T04–T15,T18–T21,T28 | 十七の厳密 QQ period 点と可逆な非 rho Jacobian |
| T22–T27,T29–T38 | 十六の正の Fp period 点と単元な非 rho Jacobian、p=11,13,17 |

QQ 証人では、完全な m 個の period 方程式に対し m 個の非 rho 変数を選ぶ。
残りを定数へ固定すると、選択 Jacobian の可逆性により実 source は rho 上 étale な
代数曲線を持つ。Poisson 接ベクトルの積分を仮定した議論ではない。

Fp 証人では全五方程式を同時に満たし、五非 rho 変数の Jacobian が単元である。
多変数 Hensel により真の Zp 解が存在し、rho=rho0+u の形式解も得られる。
従って QQ[rho] の非零多項式が admissible source 上で恒等零になることはなく、
QQ 上および C への忠実平坦な体拡大後に rho を優越する成分がある。
有限体の空性・長さや有限精度だけから標数零へ結論を外挿していない。

各証人の同じ非空開で arc 分母・固有値差・rho(rho²−1) は非零となる。
actual arc の t0 極限と対象成分の閉性から、稠密な a-image 上で Q_a が対象成分へ入る。
a↦[Q_a] の逆像も閉なので、除いた有限の極・衝突値や a=0 を含む全有限 a へ所属が延びる。
first/second rho jets と p^8 Newton 値は再現検算であり、全次数の存在は étale/Hensel 定理による。

CA、CB には別の全域 source を使える。座標順を
(h,x,y,z,w)=(old0,old1,old4,old2,old3) とし、

    f=a*h², g=−2/9,
    H=4hxw/3−hz²/3−x²z

を CA に入れる。CB では F=h²、比[a:−2/9] を使う。
両 raw tensor は divergence 零であり、全150係数で物理的 Q_a に一致する。
比は a0 を含め全域で定義される。

## 証明書と独立検証

上の四分割はそれぞれ `Q_T01_MEMBERSHIP`、`Q_MPS_EXTENDED_TORIC_ARCS`、
`Q_MPS_IMPLICIT_TORIC_ARCS`、`Q_MPS_HENSEL_TORIC_ARCS` の凍結 producer と独立 review に基づく。
最後の二件は Cloud 障害前の fresh 出力を復旧後に実読出しし、さらに全十七 QQ・十六 Fp の
完全方程式、全 Jacobian・minor、jets・Newton 更新を別の係数演算で監査した。
CA/CB の全域 lift は `Q_ALL_FAMILY_NON_CM` の独立 review に含まれる。

独立 receipt を重複なく集計した manifest は
`reports/Q_ALL_TORIC_MPS_COMPONENTS_20261103_v63.json`、実装は
`code/q_all_toric_mps_components_20261103_v63.py` である。
元の凍結証明を変更せず、全局所 primary decomposition や一般 S61 五点尽くしは未決として残す。
