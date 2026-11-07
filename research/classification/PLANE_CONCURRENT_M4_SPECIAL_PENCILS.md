# concurrent reduced cubic・m4 squarefree special j の全四 pencil

2026-11-07。主義務 G。全源有限化の数学入口であり、次元結論は別の実行証明書を要する。
rank 二 restriction を持つ plane S の zero cubic は三相異 concurrent lines、
実既約被約 D の S 沿い位数は四、最初の normal quartic は generic squarefree とする。
本ノートは I≡0 又は J≡0、かつ leading coefficient map の rank 二の場合だけを扱う。
generic j と rank 一は既受理の essential-three quintic/Saito 矛盾で空であり、再計算しない。
重根、m2/m3、非被約 R、他の rank 二支持は別残余として保持する。

## 1. 全四 leading pencil と全接続

`PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_PENCILS.md` §§2–3 の全 W 分類を使う。
I0/J0 は Witt index と singular Hankel symmetric space の common-kernel/compression により dim W≤2。
従って rank 二全残余は J0 secant、J0 compression、I0 mixed、I0 三相異の四型に尽きる。
discriminant support は R の部分集合でよく、未使用の第三線を捨てない。
base/normal の定数 GL 変換と有限 line permutation により f=xy(x−y)、
normal coordinates u,v として次の全 leading generators を取れる。

| 型 | 全 parallel leading space の generators |
|---|---|
| J0 secant | xu⁴, yv⁴ |
| J0 compression | xu³v, yuv³ |
| I0 mixed | xu⁴, yuv³ |
| I0 三相異 | (2x−y)4u³v+y(3u⁴/4+6u²v²−4v⁴) |

前三型は二つの全非零振幅を保つ。第四型は全非零共通振幅を保つ。
第四型で A=4u³v、B=3u⁴/4+6u²v²−4v⁴ とすると
I(sA+tB)=0、J(sA+tB)=4s²t−4t³、s=2x−y,t=y で J=16xy(x−y)。
したがって複素平方根を係数体へ加えず全三相異 discriminant lines を扱える。
これらへの変換は full leading pencil の normal form であり、標本点の選択ではない。

parallel Hamiltonian connection の全源性は
`PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md` の associated-graded 補題による。
squarefree binary quartic の twisted GL2 Lie stabilizer は零なので、parallel tensor が接続を一意に定める。
特に、以下に得る polynomial 接続が実際に parallel なら、off-diagonal/nonlog 接続を別に除外する必要は無い。
前三型では二つの非零振幅、第四型では共通振幅を変えても同じ接続になる。

## 2. raw P0 の完全な線形復元

f は z に依存しない。f を保つ接続の raw homogeneous representative を

\[
 P_0=\partial_z\wedge W,
 \quad W=f_y\partial_x-f_x\partial_y
       +(M(x,y)(u,v)^t)_1\partial_u+(M(x,y)(u,v)^t)_2\partial_v
\]

と取る。M は x,y 一次の全 2×2 matrix、未知は八個で、対角性は仮定しない。
各 displayed G の parallel 条件は正確に

\[
 (f_y\partial_x-f_x\partial_y)G+(Mn)\cdot\partial_nG-\operatorname{tr}(M)G=0. \tag{1}
\]

例えば y=1 chart で W の base part は 3x(x−1)∂x、normal matrix は M+f_x I。
実 divisor の Poisson identity は W(G)=div(W)G なので、normal degree 四を使うと (1) になる。
これは associated-graded parallel 条件の raw convention の直接確認でもある。
各全四型について (1) の全係数を取る八未知線形系の rank 八を保存し、一意 M を復元する。
さらに全十五個の base-linear×normal-quartic monomials に同じ operator を課し、
完全核が表の leading space に等しいことを独立 minor と全核で確認する。
これにより z を含む未検査の leading coefficient も落とさない。

この P0 の mixed brackets は z とだけ対を作り、actual D2=0 は raw {u,v}=0 を強制する。
その係数は 3fuv であり、恣意的な zero gauge の仮定ではない。
W は z に依存しないので [P0,P0]=0、P0²=0。
有限証明書は別途全150係数・全350 Jacobi・actual Pfaffian zero を保存する。

## 3. 全30/9源と最弱有限判定

既受理の normal filtration の Euler gauge により、全 tensor は P=P0+P1+P2、
P1 の全三十係数、P2 の全九係数で尽くされる。
全350 Jacobi の weight 一と actual D3=0 の全二十四行から、P1 の完全必要線形核 K1 を取る。
任意の P1∈K1 に対し weight 二の全式と全 leading-shape 条件は
L(P2)+Q(P1)=0、L は固定した同一線形 operator である。
したがって全 P1 上の fiber は空又は ker L の affine translate となる。
generic rank を exceptional fiber に流用しない。

全 compatibility quadrics を original rows の明示 left functionals と共に保存する。
必要ならその有限 rowspace の leading monomials だけを使い、真の ideal の dimension 上界を取る。
GB の完了も全 nonlinear source equality も仮定しない。
この一回の有限検査で G≤24 が得られなければ、その診断を保存し、
同じ operator の minor・jet・GB を増やさず実 D 既約性/essential-variable 幾何へ戻る。

## 4. 実在する群 fiber

concurrent 三線を保つ connected base subgroup は
z↦c z+a x+b y の三次元である（x,y の共通 scale は projective に除いた）。
normal-to-base shear 六も associated graded identity として全源を保つ。
normal subgroup は M の全同時 centralizer を使い、matrix の八 linear coefficients からその完全核を保存する。
前三型は normal diagonal torus 二、第四型は scalar torus 一である。
したがって実 connected group fiber の次元は順に 11,11,11,10。
第四型に存在しない normal diagonal 二次元を引かない。

PGL5×X→projective tensors の fiber にこの群が自由に作用するので、
固定 P0 の必要全源上界 d から像上界 24+d−g を得る。
個々の tensor の stabilizer 次元を引く議論ではない。
tensor 全体の scale は f 係数を一へ戻す projective normalization で処理する。

## 5. 証明状態

全 leading pencil 分類と接続一意性は既受理の補題に基づく数学的還元。
新しい四つの全源の次元上界は有限契約・gate 後の report を要する。
この段階では G≤24、新成分又は E 帰属を主張しない。
