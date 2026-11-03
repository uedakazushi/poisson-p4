# 元の被約 B18 交叉の conductor の独立検証

2026-11-03。著者の Lie-log source4 を再実行せず、別の有限5×5逆座標代入で
実際の N source を二変数全次数4まで作る。N19 固定空間、kernel 補償、全19
群方向、27 target 行列の左右逆を新たに照合し、全4段階の陰関数係数と全150
出力を比較する。元の25発散・350 Jacobi の全係数も独立に組み立てる。

fold の単位と distinguished 二次式は、別の5×5有理数係数問題として解く。
全出力の remainder と trace shift を新たに計算し、実際の9生成元 b_i の
phi(b_i)、phi(tau*b_i) の T³までを確認する。tau 自体に C4 値を代入しない。

## 完全な conductor の根拠

受理済みの実際の N37 座標が全完成像環を生成するため、その9個の odd
係数が I 全体を生成する。これらは abstract ν,s_i の推定ではなく、全有理
source・陰関数・Weierstrass から一意に定義される関数である。

`F=A_N×_C4 C[[T]]` の正規化は `B_N⊕C[[T]]`。B_N は A_N-module として
1,tau で生成されるので、N 側 conductor は厳密に
`J={a∈IB_N:phi(a)=phi(tau*a)=0}`、R 側は `(T⁴)` である。a∈IB_N なら
a,tau*a は A_N に属するため二条件は両方定義される。

18生成元の全 T⁰–T³ 倍からなる8×72行列と、全68次元核を独立に検査する。
QQ核だけを conductor と呼ばず、base の写像 kernel も加える。自由18加群
から IB_N への全射の前で kernel を求め、その像を取ることで完全性が従う。
像 ideal の非平坦 basechange と kernel を交換する議論はない。

## cyclic image と実際の slice

b61 の像 g、ug、tau*g、u²g は独立で全像を張り、全9個の h_i 係数も検査
する。u 作用の Hilbert 値は `(2,1,1)`。一方 B_N-module の annihilator は
7個の偶制約と
`u*tau−(116/23)u², tau²−64u²` の二次式であり、商の Hilbert 値は
`(1,2,1)`。定数 `64−(116/23)²=20400/529` が非零なので u³は従い、長さ4。
商と実像の長さが同じことにより、この annihilator は必要条件だけではない。

u は実際の中心化物理座標で `u16/2+u65−u112/2` と正確に取れる。元の全
R 曲線ではこの式が T に恒等的に一致する。19個の群 slice は偶 base の
`R_j−Rcurve_j(u)` であり、R_j は P0 で中心化した座標である。曲線を全て
含み、群微分は単位。従って群作用の形式的逆関数で N、その正規化、R、
被約交叉を同時に G19 と slice の積へ分ける。偶 base で切るため二次 cover
を保つ。群積の下で実際の像も正規化も直積になることを使い、単なる ideal
の basechange が正確だと仮定していない。

実際の `(ν,u,26物理座標)` の Jacobian 単位を照合し、19群制約を除く7座標
108,109,113,117,118,119,140 の R 曲線値が零であることを検査する。これで
上記の9生成 CI は、物理 source に束縛された完全な slice ideal となる。

`J=(b_i−h_i*b61 (i≠61), K*b61)` の証明は、任意の conductor 元を前の8差で
s*b61 に戻し、像が零となることが s∈K と同値であることによる。従って
17は**生成上界**として正しい。`J=I∩K`、最小生成数17、環の直和分裂は
主張しない。

元の D_red の全 normalization conductor をこの範囲で受理する。IB_N/J は
横断で長さ4、G19復元後に正則19変数環上 rank4自由。F/conductor は A_N/J
であり、その環全体を有限長や分裂直和としない。元の非被約 D の全環や
正確な深さには本検証を拡張しない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_Dred_conductor_review_20261103_v63.py --output work/B18-Dred-conductor-independent-review
```

60秒・2.2GB上限。著者 producer の再実行も、四次を超える新しい jet も行わない。
