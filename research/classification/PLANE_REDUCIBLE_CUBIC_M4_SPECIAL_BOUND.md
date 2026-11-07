# m4 squarefree special j：全 nonconcurrent cubic plane の G 上界

2026-11-07。主義務 G、producer `plane_collision_review`。
**結論。** rank 二 plane の zero cubic R が transverse conic+line、tangent conic+line、
nonconcurrent triangle の三型で、実 D の plane 沿い位数が四、leading binary quartic が generic squarefree なら、
\(I\equiv0\) 又は \(J\equiv0\) の全枝の PGL 飽和 projective Poisson 像はそれぞれ次元高々24,23,23。
旧非特別 j の独立受理済み上界と合わせ、三型の **全 m4 squarefree 枝** が G 非優越となる。
concurrent、重根、m2/m3、非被約 cubic、全 E は明示残余である。

幾何の全源有限化は `PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_PENCILS.md`。
有限契約 `jobs/plane_reducible_cubic_m4_special.json` の gate 後、
`code/classification_plane_reducible_cubic_m4_special.py` を一回実行した。
report は `reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_20261107.json`、
1.165687 秒で全算術 PASS。新旧 producer bytes と既存 report は上書きしていない。
gate は数学の証明ではなく、全幾何・算術の独立監査は別に必要である。

## 1. 全特殊 leading space から三つの固定 P0

前ノートでは次を全て証明した。
I0 は nondegenerate quadratic の Witt index で image rank 二以下。
その全 stable rank 二 pencil の J は三相異 root 又は二重+単純 rootで、perfect cube は無い。
J0 は全 singular symmetric space の common-kernel/compression lemma と Hankel 制約で image rank 二以下。
その全 stable rank 二 pencil は secant と compression の **両方** である。
discriminant support は R の部分集合でよく、unused components を許した。

rank 一は全て旧 scalar P0 になり、既受理の全 source 上界24,23,23をそのまま使える。
rank 二では二本以上の line components が必要なので conic+line は不可能。
I0 三相異型の三本は concurrent なので nonconcurrent triangle にも入らない。
従って triangle、\(f=xyz\) の rank 二全残余は

\[
 \langle xu^4,yv^4\rangle,\qquad
 \langle xu^3v,yuv^3\rangle,\qquad
 \langle xu^4,yuv^3\rangle
\]

という全二振幅源三つだけである。
振幅二つの非零性は generic squarefree の条件で、計算では両振幅を固定せず自由に残した。
四 root が相異なるため twisted GL2 の Lie stabilizer は零、parallel leading tensor が全接続を一意に決定する。
前ノートの diagonal connection を実際に満たすことを全五 monomial で確認し、
off-diagonal/nonlog 接続を別 parameter として捨てたものではない。

\(\tau_u=\beta_u-d\log f/3\)、\(\tau_v=\beta_v-d\log f/3\) から mixed columns を復元し、
raw normal-normal bracket は actual D2=0 で一意に復元する。
全150係数を保存した三つの P0 の raw \(\{u,v\}=k uv\) は

| 型 | k | 全350 Jacobi | actual Pfaffian |
|---|---:|---|---|
| J0 secant | 1/8 | 全係数零 | 全係数零 |
| J0 compression | 1/4 | 全係数零 | 全係数零 |
| I0 mixed | 1/6 | 全係数零 | 全係数零 |

raw \(k=0\) とする誤った仮定は使っていない。
接続の convention は affine chart z=1 でも確認出来る。
\(P_{xy}=3xy\)、\(P_{xu}=xu(1-V_{u,x}+V_{u,z})\)、
\(P_{yu}=yu(-1-V_{u,y}+V_{u,z})\) から、Hamiltonian connection は正確に
\(\beta_u=d\log f/3+\tau_u\)（v も同様）となる。

## 2. 全30/9の定数 operator による三角上界

既受理の global normal filtration により、P0 固定後の全 tensor は
\(P=P_0+P_1+P_2\)、weight 一の完全30係数、weight 二の完全9係数で尽くされる。
runner は全 basis の150係数を report に保存した。

全350 projective Jacobi の weight 一と D3=0 の全24係数を構成すると、
三型とも rank24、完全核は六次元だった。
行列、核の全30係数、独立 minor とその index を保存した。
従って全解の \(P_1\) は、正確にこの六次元線形空間に入る。

weight 二では、任意の固定 P1 に対し

\[
 L(P_2)+Q(P_1)=0
\]

という同一の inhomogeneous linear system になる。
L は全350 Jacobi と、G4 の全15係数のうち許される二 monomial 以外を零にする全13行から成る。
後者は leading symbol の全 parallel space であり、chosen point の固定ではない。
三つの \(363\times9\) operator の完全 rank は次の通り。

| 全 rank 二 leading 型 | dim P1 必要核 | rank L | 各 P1 上の P2 fiber 上界 | 固定 P0 の全源上界 |
|---|---:|---:|---:|---:|
| J0 secant | 6 | 9 | 0 | 6 |
| J0 compression | 6 | 7 | 2 | 8 |
| I0 mixed | 6 | 8 | 1 | 7 |

各任意 P1 上の fiber は空又は \(\ker L\) の affine translate なので、この上界は全 P1 に共通。
generic rank を exceptional fiber へ外挿していない。
さらに高い weight の Jacobi、actual D の既約性や squarefree の非零条件は source を増やさない。

完全な \(Q(P_1)\) と original363行の left functionals も構成して保存した。
三型とも、この段階の compatibility residual は全て零だった。
これは全高 weight Jacobi が零又は全 E 帰属という主張ではない。
今回必要な次元上界には非零 quadratic obstruction、GB、initial ideal 又は追加 jet を一切要しない。
runner に保存された quadratic rowspace は零で、その無条件な P1 次元上界六だけを使う。

## 3. 正しい群 fiber と全 PGL 像

triangle と選んだ全固定接続を保つ connected subgroup は、
base diagonal torus 二、normal diagonal 二、normal-to-base shear 六で次元10。
normal GL2 全四を引くと非 scalar 接続を動かすため誤りであり、ここでは実在する stabilizing subgroup だけを使う。
normal-to-base shear は normal filtration の associated graded で identity なので P0 を保ち、
上の完全30/9源を保つ。base/normal diagonal は自由な二振幅を保ち、
有限 triangle permutation で全 used ordered pairs を尽くす。

\(\operatorname{PGL}_5\times X_{P_0}\to\mathbb P H^0(\Lambda^2T)\) の fiber にこの10次元群が入り、
PGL5 因子への右作用が自由なので、個別 tensor の stabilizer によらず fiber dimension は少なくとも10。
固定 P0 source の全上界6,8,7から、三つの全像は

\[
 24+6-10=20,\qquad24+8-10=22,\qquad24+7-10=21.
\]

rank 一 scalar triangle の旧上界23と合わせ、triangle の全 special j 像は高々23。
conic+line は全て rank 一へ還元されたため、旧上界24,23。
これで結論が従う。

## 4. 完了範囲

I0/J0、rank 一/二、I0 の全 root multiplicity、J0 の両 stable pencil、
unused divisor components、非 log 接続を含む parallel uniqueness、全30/9源を保持した。
一回の有限 job で必要上界が得られたので、追加計算・GB・高次 jet・点列は行わない。
concurrent の失敗した旧上界31、重根、m2/m3、非被約 cubic、全 E は閉じていない。
新成分・既知成分帰属・全点分類は結論せず、G の指定全枝の非優越だけを得た。
