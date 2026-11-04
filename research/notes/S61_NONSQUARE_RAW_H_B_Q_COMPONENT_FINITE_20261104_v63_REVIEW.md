# Q成分のexact field恒等式・全22実minor・有限性の独立監査

2026-11-04。前の `CURVE_COMPONENTS_20261104_v63_REVIEW` で受理したG=Q L^4、Qの既約性、全六還元入力、Lのgeneric境界除去を用いる。今回は `COMPONENT_EXACT_QUOTIENT_Q` と `Q_ROOT_ALL_MINORS` の凍結証拠を監査する。新しいGCD、Euclid、field逆元、点、素数、GB、著者workerは実行しない。

独立native checkerの係数はQQ[e]多項式の未約分分子対(A0,A1)と分母Dであり、a²=−q1(e)a−q0(e)だけで乗法を行う。等式は分母の交差掛けで判定する。元SymPy extensionのdivision APIを使用せず、保存LC逆元はその積が1になることだけで検算する。六つのEuclid除算の全quotient/remainder、次数減少、各traceのBézout、最終三source Bézoutと全sourceの割り切れを確認する。

generic共通因子は正確に

`d+20*(900a+324e²+171e−49)/(27*(3e−2)*(6e+1)*(e+1))`

である。N=18000a+6480e²+3420e−980、D=486e³+243e²−297e−54とおけば、根はd=−N/D。保存係数とこのN/Dの等式を交差掛けで確認する。元FULLU productを同じquadratic fieldへ写した値と保存Uを比較し、U(−N/D)≠0も確認する。従ってQのgeneric候補根はLの場合と異なりFULLU内にあり、この三本だけでは除去できない。

serializedされた全係数から非定数分母を独立に収集し、保存19多項式・次数合計上界320と一致することを確認する。API失敗・修理の二つの保存receiptもhash-boundにするが、失敗を成功へ読み替えない。成功した固定証拠の全恒等式を独立演算で検算する。

## 全22元sourceの検査

実B35の同じ22 localized minor generatorを一つも省略しない。各P_iについてm=deg_d P_iとし、**Qへ還元する前に**

`E_i=D^m*P_i(a,e,−N/D)`

を独立三変数QQ辞書の積和で構成し、保存cleared_evaluationの全係数と比較する。保存されたQ quotientとa次数≤1のremainderの掛け戻しも確認する。remainder=r0(e)+a*r1(e)に対し、conjugate=r0−q1*r1−a*r1、norm=r0²−q1*r0*r1+q0*r1²を別に構成し、保存conjugate・norm・Q倍のnorm恒等式を全て検査する。非零remainderのnormも非零である。

非零indexは正確に0,2,4,6,7,8,9,11,12,13,14,16,18,19,21の十五本で、選択10,15,20を含む七本は零である。同じ検査でa^5Hのremainderは零、FULLUのremainderとそのnormは非零となる。従ってQのgeneric共通H/g/k根は全22minorの共通根ではない。

## 有限性の正確な範囲

分母例外を除いたQ上では共通根は上記の一根である。実B-bad点では十五本の非零remainderも零にならねばならず、特に一つの非零norm(e)の零点にeが制限される。このnormはeの非零多項式なので有限個である。D=0、既存field分母の零点、U normの零点も全て有限例外として保持する。

Qはaモニック2次なので、各例外eに対するaは有限個である。FULLUではa≠0、a^5Hのd4係数はa^5なので、dも有限個となる。L側の境界除去・有限分母例外と、元norm-pencilの「共通必要曲線Gと旧h/G,k/G由来の有限残余」という分解を合わせれば、**FULLU上の実B-rankbad locus全体は有限または空**である。この最後の分解は既存上流証明の入力であり、今回norm pencilを再計算したという主張ではない。

有限点を除去したとは言わない。空性、全FULLUでのB rank3、p=0でのraw H7、full cubic57、quartic239、Poisson解尽くしには昇格しない。点数をscheme lengthへ読み替えない。必要な今後の課題は保存された有限例外と旧有限残余における実rank回復・所属の判定である。

旧必要曲線と有限残余の分解・h1/k1の互いに素性は、新たに凍結された独立 `S61_NONSQUARE_RAW_H_B_ADDITIONAL_NORM_PENCIL_20261104_v63_REVIEW` を入力として結ぶ。同reviewは全norm恒等式、九GCD検算、旧1600点上界を確認している。そのgate追加は既に完了したQ/source22算術の再実行を伴わない。

専用checker、全input/evidence hash、source固定program、native logを同名REVIEWに保存する。実行バイナリはignored workに置く。
