# Actual pair gate は全 aPJ-open 上で成立する（2026-11-03）

**厳密QQの全域結果。** `a*P*J!=0` のすべての幾何学的点で、実際の二次式 Q0,Q3 を `A=x0=W=0`、B 自由に制限した2×6係数行列の rank は2である。先行 finite theorem の有限候補を、以下の保存 QQ 多項式恒等式だけで全て除外した。

これに binary-D 縮約を適用すると、補助五二次式は全 aPJ-open 上で独立となる。また、既存の gcd gate の一般点証明と全例外上の actual raw module の QQ surjectivity を組み合わせると、FULLU=`a*P*F*B*C*J!=0` 上の raw24×56 は全域で full row rank 24 である。元の product41 への17-kernel bridgeは別の監査対象であり、この packet 自体は original41/full57/ambient removal を主張しない。

自然 core 座標 `t=a/P, v=3e+20t, H=P+8t` と、既知unitによる行正規化、第2行から `h0=(160t-15H-2v+31)/36` 倍の第1行を引く操作を用いる。u²係数は `27f/40000`、`-3E/32000`、ただし `f=(v+3)(3v+4)`, `E=3H(v+4)+5v-12`。第1行の非零性、f=0 の二特殊枝の除外、三本の必要な線形t方程式と projection 恒等式は、独立 review 済みの finite theorem の全データを再検証して用いる。

f≠0 の残りで、必要な二つの elimination polynomial N,M に属する R(v) は厳密に

`R(v)=k*f(v)^2*C(v)^2*Q(v)`、k は非零QQ定数、

`C(v)=8v³+30v²+71v+58`, `deg(Q)=23`

と分解する。従ってすべての候補は C=0 または Q=0 に属する。この分解の全係数一致を検証した。Q の既約性は証明に不要である。

C=0 の場合、三本目の線形方程式 `alpha3*t+beta3=0` で `beta3=-f*C=0`。t は既知 unit なので alpha3=0。別の必要 elimination polynomial N01=0 も満たす。保存された QQ[H,v] 恒等式

`U*alpha3+V*N01=Rc(v)`

と QQ[v] 恒等式 `W*Rc=1 mod C` が同時に成り立つので、この枝に候補は存在しない。これは一つの symbolic resultant の値だけでなく、全 membership と inverse の恒等式を保存・検証している。

Q=0 の場合、保存した QQ[H,v] 恒等式

`Un*N+Um*M=L(v)*H+K(v)`

を用いる。L には QQ[v]/Q 内の保存された多項式 inverse があるので、H は保存多項式 Hbar(v)（次数22）に等しい。alpha3(Hbar,v) にも同じ quotient 内で inverse があり、t は保存多項式 tbar(v)（次数22）に等しい。この評価で三つの線形t方程式と N,M がすべて零になることも検査した。

しかし、実際の pair coefficient の第0列と u²列からなる2-minorを `(tbar,Hbar,v)` で評価すると、Qを法として保存した非零多項式 m(v) になる。さらに保存多項式 inv(v) が `m*inv=1 mod Q` を満たす。rank failure ならこの minor は零でなければならず、矛盾する。ここで用いたのは元の804項から再構成した actual pair の minor である。

この場合分けは全ての代数閉体上の characteristic-zero points に適用できる。従って pair rank failure は空である。係数環は既知 unit aPJ のみを逆にした有限型 QQ 代数なので、maximal-minor ideal の空性と Nullstellensatz、local Nakayama により corresponding cokernel も零となる。有限体からの昇格、flatness、proper compactification、追加の未知因子の逆元は証明に使わない。

再現: archive を展開し既存 Cloud activation 後に `python code/s61_nonsquare_pair_empty_QQ_20261103_v63.py`。actual804/source-unit、finite theorem の全 QQ certificate、C/Q の全場合分け、線形H membership、保存された全 inverse、最後の actual minor を再検証する。重い GB、resultant探索、factorization、CRT は再現に不要である。先行 finite theorem と bounded failure packets の bytes は変更しない。
