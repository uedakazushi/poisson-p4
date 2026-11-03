# P(1)での二CA branchの実meeting：double smooth20と全81関係

2026-11-03。個別reduced CAの二つのformal branchを対象とする。CA/CB actual intersection Dやa=0のfull ringとは区別する。

## 有限な物理curveによる二次obstruction

先行の全a≠0全旗監査により、P(1)のreduced CA completionは二つのsmooth29 branchの和集合である。
raw CA20を独立にcanonical化し、同じ物理座標における二つの150×29接行列J₁,J₂を再構成した。
literalな有理source方向u₁,u₂はcode/certificateに全29成分を保存する。
共通v=J₁u₁=J₂u₂は二接空間のintersection21に入り、実GL×a family20の接空間には入らない。

各source curveはexp(tAᵢ)(Pᵢ+t δᵢ)である。
今回Aᵢ³=0なのでexp(tAᵢ)=I+tAᵢ+t²Aᵢ²/2は**全次数で有限な実GL行列**で、逆行列とdet1も直接確認した。
Lieの二次公式を使用せず、実5×5 pushforwardの全有限座標代入を展開した。
二curveは双方とも全25 divergence係数・350 Jacobi係数がtの多項式として零であり、0次と全150の1次係数が一致する。
同一Pへ戻した2次差Lに対して

    ell=(2/5)e60+(1/10)e63+e68,
    ell(J₁)=ell(J₂)=0,  ell(L)=−5/12

を独立に確認した（全添字0起算）。従ってrank[J₁,J₂]=37にLを加えると38となる。
無限小の二次cokernel classが真に非零であり、straight方向を制限しただけの失敗ではない。
ignored workやrootのpreflight出力を実行時依存としていない。

## なぜ全meetingがdoubleになるか

実GL×a写像のrankは全a≠0で20である。
下限minor −27518828544 a¹¹/390625と、QQ(a)上の全上限rank、GL-equivarianceによってこれは全source近傍でconstant rankである。
複素解析的constant-rank theoremから、P(1)近傍の実paired-imageはsmooth20 germ Sを与える。
これは単に一つの接20平面を選ぶことではない。二つの旗liftを同時に移動するためSは両branchに入る。
代数的completionと複素解析的completionは同じ形式環なので、以下は元の代数的completed local ringについての結論である。
globalなalgebraic imageがこの一点で一つのbranchしか持たないと仮定する必要はない。

二smooth29 branchのmeeting Cの接次元は、Sの全点で正確に21である。
これはcross minor4076863488 a²⁵/3125が全a≠0で非零であり、GL作用がrankを保つことから従う。
meetingの一次式を形式的に消去すると、Cはsmooth21 ambient内にあり、Sはsmooth20 hypersurfaceとしてrho=0と書ける。
従ってmeeting ideal Kは(rho)に入り、各式はrho Gᵢと書ける。
S上の全点で接次元21だからGᵢ|_S=0。
Sが被約smoothであることにより、これは関数として恒等的に零であり、**K⊂(rho²)**である。
一点での接rankのopen性だけからこの二重因子を推論したのではない。

上の非零二次classは、ある式のrho²係数がP(1)で単元であることを意味する。
rho²×unit∈KとK⊂(rho²)を合わせてK=(rho²)。ゆえに実全meetingのcompletionは

    C ≅ C[[u1,...,u20,rho]]/(rho²).

これはP(1)およびこの証明のunit openで成立する。全a≠0への一様な二次unit主張は、この一点計算からは推論しない。

## 全CA環の81関係：kernel completeness

二branchの実imageを用いるため、全CA環はA=B₁×_C B₂、Bᵢ≅C[[29]]である。
Cの21個のcotangent generatorをそれぞれBᵢへliftし、さらに8個の一次kernel方向を選ぶと、formal inverse function theoremにより

    B₁=C[[u20,rho,y1,...,y8]],  B₂=C[[u20,rho,z1,...,z8]]

となる。各Cへの全kernelは(y1,...,y8,rho²)、(z1,...,z8,rho²)である。
このstraighteningはCが既に上の**全ring**と同定された後に行う。

共有u20を省略すると、Aの全表示は

    C[[rho,y1..y8,z1..z8,tau]] /
      ( y_i*z_j, (tau−rho²)*y_i, tau*z_j, tau*(tau−rho²) ).

関係は64+8+8+1=81個。
写像はrho↦(rho,rho)、yᵢ↦(yᵢ,0)、zⱼ↦(0,zⱼ)、tau↦(rho²,0)。
codeは全81関係の両成分での恒等式を確認する。全kernel性は次の収束するnormal formで証明する。
関係を使うと任意の形式級数は一意に

    f(u,rho)+tau*g(u,rho)+Y(u,rho,y)+Z(u,rho,z),
    Y∈(y)C[[u,rho,y]], Z∈(z)C[[u,rho,z]]

へ還元できる。tau²=rho²tauは次数を上げるので形式収束を妨げない。
両成分像は(f+rho²g+Y, f+Z)。両像が零ならy=z=0からf=0、rho²g=0よりg=0、続いてY=Z=0。
従って他のkernel関係はない。逆にCで一致する任意の二級数はこの形に書けるので写像は全射である。
これは81個のpullback恒等式だけからkernel equalityを外挿したものではない。

全81関係に一次項がないため、coneの全embedding dimensionは**38**。
以前の37は二つのbranch接空間の和であり、whole ringのedimではなかった。
0→A→B₁⊕B₂→C→0とdepth C=20より、depth A=21、dim A=29。
二つの正則minimal branchのmultiplicity加法性からe(A)=2。
射影ではsmooth factor一つを除き、dim28/edim37/depth20/e2である。

normalizationはB₁×B₂。conductorはその中で

    (y1,...,y8,rho²) × (z1,...,z8,rho²),

Aの中では(y,z,tau,rho²)であり、そのquotientは上のC、すなわち**double smooth20**である。
従ってこのgermのreduced conductor supportはSである。
この個別CAの記述をCB、actual D、a=0のfull odd idealや全global bad locusへ拡張しない。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CA_self_pair_meeting_double_review_20261103_v63.py --output <new-output-directory>

全物理curve150係数、source方向、L、cokernel functional、全81関係をgzip certificateに保存する。
既存の全旗・全rank監査を明示的に再使用する。新しいGBなし、制限120秒・2.2GB。
