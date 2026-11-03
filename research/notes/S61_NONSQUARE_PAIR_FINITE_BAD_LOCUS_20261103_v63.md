# Actual pair gate の QQ rank-failure は有限（2026-11-03）

**証明した結果。** `a*P*J!=0` 上で、実際の pair `Q0,Q3` を `A=x0=W=0`（B は自由）に制限した係数行列の rank が2未満となる集合は零次元以下である。全域で空とはまだ証明していない。従来の gcd exception の finite Nakayama 証明と binary-D 縮約を組み合わせると、FULLU=`a*P*F*B*C*J!=0` 上の actual raw24×56 の rank defect も零次元以下である。元の product41 への17-kernel bridgeはこの packet では再監査せず、元の全 product41/full57/ambient theorem は主張しない。

元の804係数項から pair chart と全12個の係数を再構成した。最初の行には保存した10項 cofactor による厳密 QQ unit identity `a*P^4` があるので、この行は open 上のどの点でも零でない。

自然な core 座標を

`t=a/P, v=3e+20t, H=P+8t`

とする。逆は `P=H-8t`, `e=(v-20t)/3`, `a=t(H-8t)`、Jacobian は P/3、`J=5H-2v-4`。従って `t,P,J` の非零性だけを使う座標変換であり、未知因子を逆にしていない。従来の QQ row normalization 後の二つの係数行をさらに P^4 で割り、第2行から `h0=(160t-15H-2v+31)/36` 倍の第1行を引く。これらは安全な row 操作である。

light monomial 順序は `(z²,v_light*z,v_light²,u*z,u*v_light,u²)`。`f=(v+3)(3v+4)`, `E=3H(v+4)+5v-12` とすると、u²係数はそれぞれ

`27f/40000`, `-3E/32000`

となる。全12係数と以下の全恒等式を certificate.json に保存した。

f=0 の二つの枝では、第1行が零でないことと rank<2 により E=0 が必要になる。従って `(v,H)=(-3,9)` または `(-4/3,7/3)`。この二つの枝の全2-minorを計算し、既知 unit t,P と非零 QQ 定数だけを除いた一変数多項式群について、保存した QQ Bezout 和が1となることを検証した。したがってどちらの枝にも pair rank failure はない。

f≠0 の chart では、u²係数と第2,3,4列（0-based）の比例条件から、t と非零 QQ 定数だけを除いた3本の必要方程式

`alpha_i(H,v)*t + beta_i(H,v)=0`（i=1,2,3）

を得る。これらが t について線形であることを全係数で検査した。二つの elimination minors から f の因子だけを除き、Hについて次数3の多項式 N,M を得た。保存された厳密 QQ 恒等式は

`C(H,v)*N(H,v)+D(H,v)*M(H,v)=R(v)`

であり、R は零でなく次数33である。さらに N を H の多項式とみた係数群に対する QQ[v] Bezout 和が1となる。従って R(v)=0 により v は有限通りであり、どの v に対しても N は H の零多項式にならないので H も有限通りである。ここでは resultant の単なる計算結果ではなく、表示した membership と content の全恒等式を再検証した。

t の係数が全て零になる場合も除外した。具体的に

`C0(v)=8v³+30v²+71v+58`, `b3=beta_2/f`, `C0=-beta_3/f`

とおくと、保存した QQ[H,v] 多項式によって

`U*C0+V*b3+Z*alpha_3=1`

となる。この証明には係数が `beta_2/f` に現れる `20v²-20v-115` を点ごとに逆にする操作はなく、一変数 Euclid から導いた多項式恒等式のみを用いる。仮に全 alpha_i が零なら必要方程式から全 beta_i も零となり、上の和1に矛盾する。従って各候補 (H,v) 上で少なくとも一つの alpha_i が非零であり、t は高々一つ。これで pair failure の零次元性が従う。

この証明は finite QQ theorem であり、32003での unit 計算を QQ へ昇格させていない。全 pair failure の空性を得るための proper finite Z_(p) model、全 localizer の unit 性、または literal QQ unit identity は未構成である。

再現は既存 Cloud activation 後 `python code/s61_nonsquare_pair_finite_bad_locus_20261103_v63.py`。約6秒で actual804→pair係数、QQ source unit、自然座標、全線形t方程式、resultant membership、H-content 和1、t-coefficient 和1、両特殊枝の全 QQ 恒等式を検証する。重い GB・追加有限体・CRT を用いない。
