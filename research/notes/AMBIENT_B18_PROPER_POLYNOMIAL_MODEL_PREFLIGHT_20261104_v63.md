# 実 E4 の proper-source 多項式模型: 14変数・17式の定数 QQ preflight

2026-11-04。これは normalized-source tensor E4 の模型であり、image-kernel で定義される
実 D0 と同一視しない。新 source jet、GB、多項式 parameter 行列の展開は行わない。

## 1. affine/scalar と固定 source の範囲

非零の生 tensor P0 を扱う affine cone の局所模型である。
projective P18 fiber にする追加 quotient は取らない。
P0 の全 proper CA/CB fiber は一つの幾何点を持つ長さ2であり、その近傍では
smooth proper source が実 normalization であることが既独立受理済み。
今回の source はその同じ閉旗の局所 chart を使う。

零 weight torus の fixed locus では定 vector line L=ey、CA の g=0、CB の ratio mu=0。
全 source 上の L 移動四方向と g/mu は正 weight で、fixed locus に残さない。
残る旗は covector ell と二次元 W の旗 F(1,2;4)、五次元である。
各固定旗の生線形空間は f三方向+H十六方向、十九次元。
従って両 normalized zero source はそれぞれ5+19=24次元。
全 proper fiber と source chart の受理済み有限性を用い、Bruhat の dimension count だけから推定しない。

物理座標 (h,x,y,z,w) に対し reference旗は ell=x,W=span(x,h),L=ey。
global source の (markedh,k,y,z,w) を (x,h,y,z,w) に交換し、全150 canonical 係数を独立構成した。
CA は raw(f,0,H)、CB は rawB(f,0,H) の十九基底。
中心の f=-markedh²、H=-markedh*k*w+k²*z+z³/3 を元十九係数へ戻し、
二つの full150 点が既受理 physical P0 に厳密に一致することを確認した。

## 2. 対角 GL4 の十五方向と二種類の実切断

GL4 は (h,x,z,w) に作用し y を固定する。
GL4 の中心はこの二次 ey-wedge tensor に -P0 の微分を与えるため、tensor scalar は既に含まれる。
y-scaling と別の projectivization を数えて次元を一つ減らしてはいけない。
全16 Lie方向の物理 rank は15であり、全GL4ではなく十五変数の formal group section を使う。

reference旗を保つ parabolic GL4 は、row x が column x だけ、row h が columns h,x だけ、
rows z,w が自由という十一 Lie方向を持つ。
その CA 十九 fiber への実係数行列は rank10。
独立な十方向と十 fiber 座標を QQ で選び、残る九 fiber 座標を自由とする。
これで sourceCA旗五変数を零にし、さらに十 fiber 係数を中心値へ固定する。

残りの五 Lie方向は E_xh,E_xz,E_xw,E_hz,E_hw。
source旗座標の微分はそれらに対して -I5。
その旗を reference frame へ戻す操作が同じ bad Lie方向を打ち消すため、
選択十 fiber 座標への微分は零となる。
parabolic十方向に対する fiber 応答は上で選んだ10x10 unit 行列。
従って source gauge15 の実ヤコビ行列は block diagonal(-I5,M10)。

さらに同じ十五 GL4 方向を**既存の物理 common zero cut15**へ代入した全15x15行列も
QQ 上で構成し、det が非零であることを確認した。
二つの具体的 unit 行列と全 source 基底は添付に保存する。
この物理 cut は current E4 の零 weight 切断そのものであり、四つの positive cut は
fixed locus 上で自動的に零になる。

## 3. all-order gauge 同値の証明

二つの proper zero sources の積は smooth48次元で、物理 tensor 等号の incidence を持つ。
上の source gauge15 は一方の CA source だけを使う十五個の正則条件。
十五変数 GL4 section と source gauge slice33 の action map は、上の unit15 により
全 source48 germへの形式同型である。
incidence は GL4 不変なので、その pullback は group15 と slice incidence の積となる。

その積の中で、既存 physical cut15 は group15 に関して別の実 unit 行列を持つ。
従ってこれらを unit IFT で解くことで、source-gauge incidence と physical-cut incidence は
全 scheme として形式同型となる。正則解を高次で推測していない。
残る stabilizer 一次元を quotient していないし、reduced support や tangent rank の一致だけを
使った同一視でもない。
群 section は elementary (I+tEij)、対角 (1+t) の積として代数的に取れ、全分母は原点の単位。

## 4. 十七の literal 多項式条件

source gauge 後の CA fiber は中心+九自由基底の affine slice。
CB の相対旗を

    ell=x+a*h+b*z+c*w,   k=h+d*z+e*w

とし、物理 reference順に並べた coordinate frame を

    hbar=k, xbar=ell, zbar=z, wbar=w, ybar=y

とする。determinant は正確に1。
逆は h=hbar-d*z-e*w,
x=xbar-a*hbar+(a*d-b)*z+(a*e-c)*w なので parameter degree は高々2。

元の canonical zero36座標における固定 CB十九空間を定数行列 B とする。
その十九 unit rows と左逆を選び、完全な17x36 annihilator N_B を構成した。
N_B B=0、rank N_B=17、kernelは固定CB十九空間そのもの。
CB fiber の十九未知係数をこの左逆で消すと、残る incidence は正確に

    N_B * coordinates36(frame_* (P0 + sum CAfree9_i*v_i)) = 0.

これは十四変数 (CAfiber九個+相対旗五個)、十七多項式の**全 ideal**表示である。
係数も frame も全て添付に固定され、未指定の関数はない。
この packet ではその170個の parameter coefficient 多項式をまだ展開していない。
CBsource係数を消す左逆は定数 QQ unit のみで、conductor・p・g を反転しない。

五つの relative flag derivative と九つの fiber derivative から得る実17x14ヤコビ行列も構成した。
その rank は5、非零5 minor と pivot rows/columns を添付に固定する。
従って五 unit を消した embedding dimension は9で、current E4 と実 source gauge で対応する。
これだけで残り ideal の CM、graph 冗長性、socle は結論しない。

## 5. 一回 builder の budget と停止条件

ey を固定する frame 作用は四変数 vector field の一つの contravariant factor と
二つの inverse coordinate factors だけを持つ。五旗 parameter の次数は高々1+2+2=5。
CAfiber係数へは線形なので、全十四変数での次数は高々6。
十個の CA列 (中心一列+自由九列) に対し、17x10個の係数多項式があり、
各々五変数次数五以下なので最大252項、全42840項。
これは150000項・32MiB上限以下である。
定数分母の lcm と保守的な分子 bit 上界も report に固定する。

次に許可を得た一回の builder は、この全170多項式を展開して全CB左逆・元36座標との
恒等式と上の constant Jacobian を照合するだけでよい。
60秒・2GB・32MiB・150000 expanded terms を超えたら入力と停止状態を保存し、
GBや高次 jet を自動追加しない。
current四 core への座標変換や残る十ニ式の syzygy は、その後の別問題である。

再現入口:

    .venv/bin/python code/ambient_B18_proper_polynomial_model_preflight_20261104_v63.py --output work/B18-proper-source-preflight

この入口は定数 QQ の preflight のみで、170 parameter 多項式の展開はしない。
