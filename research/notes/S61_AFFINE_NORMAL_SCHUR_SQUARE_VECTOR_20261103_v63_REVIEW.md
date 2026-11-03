# controlled Schurのvector／square条件と次数staircase：独立検証

`S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63` の二つの正確な空性条件を受理する。併せて `S61_AFFINE_NORMAL_SCHUR_STAIRCASE_20261103_v63` の実際の3350 coefficient-entryと普遍companionの次数証明を独立に確認し、formal45変数でのtarget次数を7／14へ改善する。**実際のmembership、unit、一般五点尽くしは未計算**である。

## 固定する実algebraと仮定

generic betaの係数体と新しい独立moment変数aを含む体をkとする。formal coefficient環は

\[
R'=k[n_0,\ldots,n_{41},S_1,S_2,H],\qquad
\mathfrak m=(n_0,\ldots,n_{41},S_1,S_2,H)
\]

である。S1,S2,Hの実Aへの作用は既知のnormal polynomial circuitで決まり、自由な別source parameterではない。全134式を保つglobal quinticから

\[
R'^{670}\xrightarrow{\Phi}R'^5\longrightarrow A\longrightarrow0
\]

という実presentationを使い、target基底は **1,h,h²,h³,h⁴** である。したがってe0は実unitに写る。

既存の独立受理済みsharp moment proofにより、四つのmoment columnsは\(\operatorname{im}\Phi\)内にあり、pivot \(\delta\) はknown E5上で零、未知factor B(a)上で非零因子である。associated componentsやnilpotentsを除いていない。kへのfaithfulな係数体拡大も確認済みである。このreviewはそのactual source／norm proofをhash固定して再利用する。

## 正確なkernelとvector条件

仮にB=0ならA=E5である。hの五つの冪はk-basisなので、

\[
\ker(R'^5\to E5)=\mathfrak mR'^5
\]

が成分ごとに厳密に成り立つ。E5が五つのk-rational点へ分裂する仮定は不要である。

\(\delta\in\mathfrak m\) だから、この場合 \(\delta e_0\in\operatorname{im}\Phi\) である。逆にこのvector membershipが成り立つと、Aで \(\delta\cdot1_A=0\) となり、特に \(\delta\cdot1_B=0\) である。B(a)での非零因子性からB(a)=0、faithfulnessからB=0となる。よって

\[
B=0\quad\Longleftrightarrow\quad\delta e_0\in\operatorname{im}\Phi.
\]

\(\delta\)をA全体の非零因子と仮定してはいない。known factorを殺すため、それは通常零因子である。

## scalar square条件

四pivot列を \(V=\binom{u}{D}\) とし、

\[
\lambda=(\delta,-u\operatorname{adj}D),\quad\delta=\det D,
\qquad J_g=\lambda(\operatorname{im}\Phi)
\]

とする。B=0なら上の正確なkernelにより

\[
J_g=\mathfrak m(\lambda_0,\ldots,\lambda_4).
\]

\(\lambda_0=\delta\in\mathfrak m\) だから、\(\delta^2\in J_g\) となる。

逆に、各 \(\lambda\Phi_j\) は四pivot列と元の一列による5-minorである。pivot列が元Φの係数体上の定数結合であることとCauchy–Binetから、\(J_g\subset\operatorname{Fitt}_0(A)\subset\operatorname{Ann}_{R'}A\) である。従って \(\delta^2\in J_g\) ならB(a)をannihilateし、その非零因子性からB=0となる。よって

\[
B=0\quad\Longleftrightarrow\quad\delta^2\in J_g.
\]

これは\(R'\)上のglobal membershipであり、δを反転した環だけのmembershipではない。旧cubic-power条件も有効である。これらの存在証明からcofactor次数の上界を与えない。

## 独立なtoy検査

20独立変数のgeneric5×4行列に対し、別実装のsymbolic determinantで五つのcofactor、全四つの \(\lambda V=0\) を確認した。100列 \(x_j e_i\) のimageは \(\mathfrak mR'^5\) であり、\(\frac14\partial\delta/\partial x_j\) を明示的に用いてvector \(\delta e_0\) とscalar \(\delta^2\) の両membershipをQQ polynomial identityとして検査した。これはuniversal module proofの検査であって、実S61 targetのmembership標本ではない。

## 新しい次数7／14／30の独立なbinding

保存された全134 remainderから、producerのmonic divisionとは異なる一段のcompanion recurrence

\[
(c_0,\ldots,c_4)\mapsto
((H-g_0)c_4,c_0-g_1c_4,\ldots,c_3-g_4c_4)
\]

を用いて、全670列×5 entryをQQ疎多項式として再構成した。全3350 entryについてdegree≤2、初期h1,...,h4係数はlinear、quadratic部分はtarget row k>jで零、k=jで元のQ_iに一致するという全staircaseを検査した。実際の非零Qと非零r4のwitnessも保存値に一致する。

四moment列j=0,...,3のbottom h4行は全てlinearである。他の三行の次数は2以下なので、\(\deg\delta\le7\) である。普遍companionでは、degree7を取り出す際にh1,h2,h3行をquadratic部分、h4行をlinear部分へ置き換える。そのdeterminantを独立に計算すると厳密に \(-L_4Q^3\) であり、sourceの二つの非零witnessからgeneric degree7も確認できる。実際の巨大δは展開していない。

Cramer numerator Uはh1行をtop行へ置き換えたminorで、依然linearなh4行を含むためdegree≤7である。全cofactor λ_iの次数が7以下だという主張ではない。

従ってformal45 parameter次数に関して、vector targetは7以下、square targetは14以下、旧cubic targetは21以下となる。135式のalternativeでは、\(\delta^4 r_i(U/\delta)\) の係数次数2と四つの次数7因子から、134式の次数は30以下であり、norm式は以前どおり10以下である。

これらは非線形のS1(n),S2(n),H(n)への代入後のnormal次数ではない。cofactor degree boundや実membershipを意味しない。旧8／16／34という緩い上界もそのまま有効である。

独立codeは `code/s61_affine_schur_square_review_20261103_v63.py`、receiptは `reports/S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63_REVIEW.json` である。両producerの新規再実行で証拠のbyte一致も確認する。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/s61_affine_schur_square_review_20261103_v63.py \
  --work work/Schur-square-staircase-independent-replay \
  --output work/Schur-square-staircase-independent-review.json
```
