# reduced reducible cubic plane：全接続入口と parallel normal symbol

2026-11-07、主義務 G。E優先作業への切替時点の確定数学を保存する。
**G次元上界は未完了。** R非被約は別補集合。四被約型を単に log と仮定していない。

## 全有限接続源

四型を `z(xy−z²)`, `z(xz−y²)`, `xyz`, `xy(x−y)` とする。
H0(Ω1_P2(3)) は八次元。tracefree linear field V と三次元 volume vol に対し
`αV=iE iV vol/f` がその全源である。
tracefree field三つ H,A,B に対応する全24係数接続の flatness は正確に

```
H(f)+det(E,A,B)=0,
A(f)+2det(E,H,A)=0,
B(f)−2det(E,H,B)=0.
```

各三次係数を取った全30式であり、全三成分の曲率を別途展開して一致させた。
`CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_CONNECTION_20261107.json` は全24源・全30式・
全log/closed源と、全四型それぞれで真正 nonlog flat 接続を保存する。
有限 gate 後0.0632秒、RSS13,824KiB、全恒等式 PASS。

logarithmic 部分だけを扱うと、全log形式次元は順に1,2,2,3、closed scalar形式次元は1,1,2,2。
flat sl2部分の次元は3,3,4,6。
接線conicでは basis α,β を `dα=0,dβ=η,α∧β=2η` と取れ、
`Ω=αT+βU` の全条件は `2[T,U]+U=0`。
triangleでは二行列が commute、concurrentでは二つのclosed formの行列は任意で第三係数だけが零。
これらは log 部分の結論であり、全接続とは主張しない。

真正 nonlog witness は線形2×2行列 M から
`Ω=M⁻¹dM−(1/2)dlog(detM) I` として作った。
conic二型は symmetric M、三直線二型は `M=[[x,z],[0,y]]`。
全曲率零・horizontal・simple pole と、あるcomponentの非log polar partを exact QQ で確認した。
従って irreducible cubic の π1=Z3 による log 化はここでは使えない。

## 最初の normal symbol の全場合に成り立つ補題

S=P2⊂SingD、π|Sのrankが二、zero cubic f が reduced とする。
S沿いのDの位数mは2,3,4。m=5ならDは二normal変数だけのbinaryquinticとなり複素数体上で分解する。
最初の項を
`Gm∈H0(O(5−m))⊗Sym^m C²` と書く。
adjunction `detTP4|S=detTS⊗detN` とπSで割って得る

`γm=Gm/f ∈ Sym^m N*⊗detN`

はU=S\R上で conormal flat connection の parallel tensor である。
証明は Hamiltonian flow がπ²とπSを保存することをI-adic associated gradedへ制限するだけでよい。
{I,I}⊂I²によりnormal上の作用はbase functionのliftによらず、Jacobiからflat。
Hamiltonian fieldsがTUを張るためγmは全方向にparallelである。
log性・ordinary-double仮定・formal splittingは要らない。

generic binary formがsquarefreeなら、そのdiscriminantはUで零にならない。
従ってΔ(Gm)の零divisorはRの成分だけに支えられ、

`Δ(Gm)=const*product(fi^ai),  sum(deg(fi)*ai)=(5−m)(2m−2)`。

右の次数はm=2,3,4について6,8,6。指数は有限で、全非log枝も含む必要条件である。
m=2ではtwisted GL2作用でdiscriminantが不変なため、更に `Δ(G2)=const*f²`。
m≥3でrootが三つ以上なら相対作用のLie stabilizerは零なので、connectionはγmから一意に有理復元出来る。
generic repeated-root partitionは別の明示補集合として保持する。

受理済み型7第五23成分は、plane `{h=x=0}` で
`R=10*y*z*(y−81z)`、m=3。
leading binary cubicのdiscriminantは `−12*y³*z³*(y−81z)²` で上の有限指数条件と一致する。
したがって三concurrent型の全不存在は偽であり、弱いG次元上界が目標である。

## m4での追加有限化と有限modeの結果

binaryquarticのclassical invariant I,J が両方非零なら、parallel性よりJ²/I³はconstant。
G4係数はbase一次なので、UFDから `I(G4)=aL²,J(G4)=bL³`、LはRの一次因子。
squarefreeならΔ(G4)もL⁶であり、G4のprojectiveconnectionはL以外では正則。
元のsimple-pole仮定とH0Ω1(L)=0からtracefree接続は零。
よってG4=L*constant binaryquartic、scalar接続差は `−(1/6)dlog(L³/f)`。
I=0/J=0、非squarefree及びm2/m3は残る。

reusable runner `classification_plane_reducible_cubic_connection.py` の
`--scalar-fourth-order` は有限契約後、この四P0に全30weight1/全9weight2を課した。
結果V2は `CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SCALAR_V2_20261107.json`。
位数四条件D3=0を含むweight1kernelは8,8,8,14、weight2fiber次元上界は3,3,3,6。
従って最初の三型の全G上界は24,23,23となり、このm4squarefree・I,J非零枝は閉じる。
concurrentは粗上界31で未完了。
二回目にはG4=L*g4という全10leading-shape条件と完全なP1二次polarizationも加えたが、
P2線形rankは改善しなかった。この場合の真の制限はP1二次compatibility側であり、
単にlinear fiber rankを増やす経路は二回で停止した。GB/jet farmへは進めていない。
詳しい矢印と補集合は `PLANE_REDUCIBLE_CUBIC_M4_SCALAR_BOUND.md`。
default初版とm4初版の両reportはembedded sourceを含め不変保存する。
