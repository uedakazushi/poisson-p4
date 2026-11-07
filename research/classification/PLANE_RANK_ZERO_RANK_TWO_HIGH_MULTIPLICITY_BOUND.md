# rank 零 plane・rank M 二・高重複度：閉じた全源と三つの残余

2026-11-07。主義務 G、producer `plane_collision_review`。
**到達点。** 実既約被約 quintic D の rank 零 plane、k=1、rank M 二、m≥3 の全 first-symbol 源を、
`PLANE_RANK_ZERO_RANK_TWO_HIGH_MULTIPLICITY_FIRST_JET.md` で全補集合付きで還元した。
M(p)非零の全四源は projective 像次元≤23。
M(p)零の全源は、ND の全 row と Borel の大部分が rank≤2 のため空、
Borel の一つの wall は actual Pfaffian の二次因子のため空となる。
**残るのは Borel \(r=(\rho,0)\)、\(\rho=1,5/3,3\) の三つの固定 first symbol だけ**。
これらは依然未解決であり、全 high-m G の完了とはしない。

## 1. 一回の全30/9計算と不変証明書

新実行器 `code/classification_plane_rank_zero_rank_two_high_m.py` は shared runner を import・編集せず、
既存 Fraction spatial engine と新しい exact univariate linear algebra を使用する。
契約 `jobs/plane_rank_zero_rank_two_high_m.json` の gate 後、4.220093 秒、全算術 PASS。
不変 report は `reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_HIGH_M_20261107.json`。
その `all_G24_closed` は **false のまま** であり、粗 linear bound だけで未解決源を閉じていない。

全七新モデルの全150 P0係数、全350 Jacobi/Pfaffian 零、全350×30 D1、全350×9 D2、
完全 kernel、nonzero minor、全 univariate factor/root を保存した。
全39 higher coefficients を保持し、任意 P1 上の P2 解は同一 D2 の inhomogeneous fiber なので
\(d\le\dim\ker D_1+\dim\ker D_2\)。
NDSS0 は旧 nonkernel report の NDSS s0 を新しい幾何で復活させた。
同 root の全行・kernel/minor は旧別作者 review でも既に検証されているので、新たに再計算していない。
旧 m二 theorem の scope 自体を広げたのではなく、その不変の算術を新しい source 上界へ適用している。

## 2. M(p)非零の全四源

| 全 high-m source | affine symbol dimension | generic d | exceptional d 最大 | 全像次元上界 |
|---|---:|---:|---:|---:|
| NDSS0 | 10 | 8 | — | 23 |
| NDN0 | 9 | 8 | — | 22 |
| BSS0(τ) | 10 | 7 | 8 | 22 |
| BN0(τ) | 9 | 7 | 10 | 23 |

BSS0 の全 selected-minor roots は
\(0,-3/2,1/6,1/2,3/2,-1/2\)。
τ0,1/2 で d8、他で d7。
BN0 の全 roots は \(1,0,-3/2,1/2,-1/2,3/2\)。
d は順に9,8,9,10,8,10。
全 root は QQ 上の exact certificate を持つ。
零 parameter も明示的に含め、generic open と全 root が C 全体を覆う。
fixed parameter では symbol base が一つ下がることと plane choice6・projectivization−1 を戻した。
したがって全四源が弱い G 非優越を満たす。

## 3. R0 の巨大な linear kernel は実際には common-U 型

初回の粗 d は18を基本とし、像上界31–33以上になった。
しかし完全 kernel の係数を読むと、次の全源ではより強い結論を得る。

* ND 非零 row \((1,\rho)\) の全ρ（selected minor の唯一の root ρ0も含む）。
* ND row零。
* Borel rc非零の全 open source。
* Borel rowline の generic minor open と roots \(0,-3,-1/3,1/3\)。

いずれも \(\ker D_1\) は十二次元で、その全 tensor は

\[
 \Pi_1=\partial_a\wedge V_1(b,c,z,w)
\]

と書け、係数は a に依らない。
全30源の非零 row index は正確に
\(1,2,4,5,10,11,13,14,18,21,24,27\) だけである。
report の basis order は旧30源と同じであり、この index は全 kernel の各 column で確認出来る。
\(\ker D_2\) は六次元で、全て
\(\partial_a\wedge\partial_b\)、\(\partial_a\wedge\partial_c\) に normal quadratic を掛けたもの。

\(\Pi_0\) も \(\partial_a\wedge V_0\) で a-independent なので
\([\Pi_1,\Pi_1]=0\)。従って weight 二 equation の forcing は零で、
実際の全 \(\Pi_2\) は \(\ker D_2\) に入る。
全 tensor が \(\partial_a\wedge(V_0+V_1+V_2)\) となり、projective rank は高々二。
rank 四 actual D の枝としては空である。
粗 tangent dimension を改善したのではなく、全 higher solution の common factor を証明した。

## 4. Borel rowline ρ=−1 の actual 二次因子

Borel rowline の残る exceptional parameters は \(-1,1,5/3,3\)。
ρ=−1 では D1 kernel は十三次元だが、\(\partial_a\) を含まない bivector direction は
\(\partial_b\wedge\partial_c\) だけである（全 source row7,17）。
全九 \(\Pi_2\) も base pairs だけなので、この場合は nonlinear equation を解く前から

\[
 \Pi=\partial_a\wedge V+F_2\partial_b\wedge\partial_c
\]

となる。\(F_2\) は homogeneous quadratic、任意の higher coefficient を含めてよい。
従って \(E\wedge\Pi^2\) には F2 が因子として入る。
F2零なら actual Pfaffian零、非零なら二次×三次であり、既約被約 quintic に矛盾する。
\(\Pi_2\) を linear kernel と仮定する必要も無い全 source 因子証明である。

## 5. 最後の三つの全 source と一回の flag criterion

ρ=1,5/3,3 で D1 kernel はそれぞれ14,14,15、D2 kernel は六。
fixed-symbol coarse image upper はそれぞれ32,32,33にすぎず、ここを G≤24 とは扱わない。
全 source を保ち、Borel の自然な normal flag \(w=0\) に対して一回だけ次を調べた。

\(\Pi_1\) の完全核の全 parameter と \(\Pi_2\) 全九係数を残す。
全350 projective Jacobi の parameter 次数≤二の係数行空間へ、
actual Pfaffian の modw 全 quintic 係数が入るかを QQ で判定する。
成功なら original Jacobi の定数 linear combination が literal な Pfaffian factor を与える。
これは一般の ideal membership や高次 multiplier を含む GB 問題にはしていない。

専用契約 `jobs/plane_rank_zero_rank_two_high_m_flag.json` を gate し、
同 reusable runner に別 mode を加えて1.127933秒で実行した。
旧十 FunctionDef AST は初回 report の embedded source と完全一致、
旧 report/producer bytes は不変で、初回 job の runner binding だけ理由付きで refresh した。
新証明書 `reports/CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_HIGH_M_FLAG_20261107.json` は
`all_three_whole_flag_exclusions_pass=false` を保持する。

三型とも **唯一 z5係数だけが定数 Jacobi rowspace に入らなかった**。
それ以外の全 modw 係数は、保存した exact original350-row functional により零となる。
従って任意の実 Poisson solution に対し

\[
 D|_{w=0}=\kappa z^5
\]

という全 parameter の幾何条件は証明されたが、κ零は未証明。
非零残差はρ1,5/3で parameter monomial (8,20)、ρ3で (8,21) の負号付き一項であり、
全 monomial labels、行列、成功 functional、残差を保存した。
一方の source を点に特殊化した結果ではない。

二つの bounded criterion はここで終了する。
同じ route の multiplier degree、GB、jet、点列を増やさず、
以後はこの三全 source の新しい幾何・明示 family 又は適切な次元証明が必要である。

## 6. 完了と残余の正確な区分

M(p)非零は全四源で G≤23。
M(p)零の ND 全 row、Borel rc非零、Borel rowline の三点を除く全 parameter は空。
これには零 row、nilpotent axis、全 minor exceptional roots、ρ−1の因子型を含む。
残る三点は必要条件 \(D|_{w=0}=\kappa z^5\) を持つが未解決。
初回/二回目 report の false を消さず、閉じた部分と非閉鎖を分けて保存した。
rankM三、rankM≤一、k0、全 E と global/main/Git は本作業で変更しない。
