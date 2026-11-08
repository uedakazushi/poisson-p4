# nil31短鎖 r=2/3：二つの canonical SS12 companion placement の全源障害

2026-11-08、PRIMARY E、`/root/quartic_review` の純証明。
短い Jordan block を SS12 の小さい/大きい weight に置く二つの明示 companion を扱う。
一般の SS12 modular deformation がこれらの regular normalform になるとは証明していない。
従って任意 SS12 閉包からの非所属・wholeE完成を主張しない。
第一 companion は gate 後一回の exact fullsource diagnostic、第二は新計算無しの手証明である。
新 profile、matrix、jet、GB の実行は追加しない。

## 1. 全 target と全 SS12 parent

target は canonical(15)、t!=0、d任意。
特に旧座標の Pxy coefficients は hy に −t、zw に2t、z²/yz に零である。
実source全六parameters/正規化inverseと閉包の補集合は元支持ノートのまま保持する。

SS12 source は全90線形方程式の九parameter解であり、全 nonlinear Poisson 条件は

    g1=n4(n7+n8−1),g2=n2n4,g3=n1n4−n3(n7+n8),
    g4=n0(n7−1),g5=n6(n5−n8),g6=n2(2n5−1),
    g7=n0n6+n1(1−n5)−n2n3。

全九を free Poisson space と呼ばず、特に g1/g4 の product relations を保持する。
以下 n_i は任意 meromorphic functions、有限 ramification も許す。
finite Taylor ansatz に限定した排除ではない。

## 2. 短 block weight一の companion

modular path と全 spectralframe は先行 diagnostic に固定した。

    Nε=y∂x+z∂y+(-4ε³x+4ε²y−εz)∂z+(1+εw)∂w,
    X=z−εy−2ε²x,Y=z−4ε²x,Z=h+εw,W=z+3εy+2ε²x。

全sourceは εFε^-1_*P_SS(n)。旧 Pxy の hy、Pxz の hz が全九方向＋baseで恒等零。
target は −t,+t なので矛盾する。
全 negativepowers/150endpointの298×63線形系には、unit row functional二本が保存されている。
任意 tails の排除は先行手証明の monomial欠落そのものからも従う。

## 3. 短 block weight二の companion と actual inverse

もう一つの canonical spectral placement は

    Nε=y∂x+z∂y+(-2ε³x+ε²y−2εz)∂z+(1+2εw)∂w。

long block eigenvalues は(-2,-1,1)ε、short nonzero eigenvalue は2ε。
hを固定した全 spectralframe は

    X=z−ε²x,Y=z+εy−2ε²x,Z=z+3εy+2ε²x,W=h+2εw。

inverse は

    x=(2X−3Y+Z)/(6ε²),
    y=(-4X+3Y+Z)/(6ε),
    z=(8X−3Y+Z)/6,w=(W−h)/(2ε)。

generic ε!=0 で可逆、全tensorに同じ pullback を施す。
任意 regular endpoint は全 SS12 source の εFε^-1_*P_SS(n) として読む。
ε0で inverse が可逆とは仮定しない。

long derivatives の xy wedge 成分は正確に

    (∂X∧∂Y)xy=−1/(6ε³),
    (∂X∧∂Z)xy=+1/(6ε³),
    (∂Y∧∂Z)xy=−1/(6ε³)。

W derivative は旧 ∂w のみなので、W-index を含む括弧は旧 Pxy へ寄与しない。
long three brackets を全て戻すと

    PXY=(n5+n7−1)XY,
    PXZ=n0hY+(-n5+n8+1)XZ,
    PYZ=n3h²+n4WX+n5YZ。

したがって旧 Pxy の coefficient は厳密に

    [hy]Pxy=n0/(6ε),
    [zw]Pxy=−n4/(3ε)。                              (1)

他の全 parent parameters はこれらの monomialへ寄与しない。
target は −t,2t なので

    n0/ε→−6t,n4/ε→−6t。

t!=0 により n0/n4 は generic field で非零。
**元の** product equations g1/g4（加えてg2）から

    n7=1,n8=0,n2=0                                  (2)

が恒等的に従う。generic period を独立Taylor係数へ自由化したものではない。

## 4. 相反する二つの regularity residues

次に同じ旧 Pxy の z²/yz coefficients を取る。
Y/X/Z は全 long variables、WX は h 又はwを含むのでこれらへ寄与しない。
直接 product expansion から、(2)の下で

    [z²]Pxy=(1−3n5)/(6ε²),
    [yz]Pxy=(3−8n5)/(6ε)。                           (3)

例えば XY の z²/εyz coefficients は1/1、XZ は1/3、YZ は1/4で、
上の三 wedge signs と全 n5 coefficientsを戻した式である。
他の source coefficientsを零に選んだものではない。

regular endpoint だけでも第一式は n5→1/3、第二式は n5→3/8 を強制する。
同一 meromorphic function に両方は不可能。
targetの二係数は零なので更に強い消失次数も必要だが、矛盾には用いない。
d、磁気定数、higher tails、全他sourceを保持したまま、t!=0全targetを排除した。

## 5. 使える結論と保留する範囲

両 canonical companion path について、全 SS12 quadratic source の任意 parametercurve は
shortnil r2/3 target に regularに到達しない。
第二の排除には新数値実行・新operator・追加Taylor次数を一切用いていない。
第一の nonlinear support計算は実行前に止まり、第二は true nonlinear products 二つと
元coefficientsだけで閉じた。

任意の SS12 spectralflag deformationをこの二companionへregularに移せる定理は無い。
coupled long/short flags、変動 eigenratio、異なる modular type、repeatedHを失う親全族は残す。
正の actualcurve又は別の knowncomponent recognitionが全r2/3 E課題の次の義務である。
同じ二companionの係数次数・profiles・jetsは増やさない。
