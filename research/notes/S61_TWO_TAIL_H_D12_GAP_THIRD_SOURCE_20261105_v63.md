# H上の二tail十分条件のgapと固定第三tail：条件付きsource

2026-11-05。D12商環の非零一次結果は独立値peerを待つ。本noteは条件付きsourceのみ。
新しいtailや行列式の値はまだ計算しない。

Hは既受理のsquarefreeな厳密共通GCDで、完全な反転積U=adFBCJと互いに素である。
D12の分母はUの因子だけなので、D12=N/uと書く。ただしuはUの積の単位である。
独立peerによってHがD12を割らないことを受理したと仮定する。
もし複素幾何点でV(H)∩D(U*N)が空なら、Nullstellensatzとchar0 geometric squarefreeにより
U*Nはsqrt((H))=(H)に属する。UFDでgcd(H,U)=1を消すとH|N、したがってH|D12となり矛盾。
よってV(H)∩FULLU∩D(D12)は非空で、Hの少なくとも一つの成分の非空な二次元openを含む。
これは選択二tailがともに消える場所であって、実bad Fittingやrank低下の主張ではない。

## 既存sourceから固定する第三候補

既存TAIL8一覧で、未計算のoriginal32 target3を一つだけ選ぶ。
target1とtarget2は既に計算済み。残りではtarget3だけが同じ149-node列閉包を持ち、
target4/6/7/8/9の172-nodeより小さい。randomなentry/minor探索ではない。
pair行はcurrent15 `[2,4,5,8]`、旧Q2 rootsは `[2551,2761,2866,3181]`。
既存cached kappa四個は旧root `[179,198,217,237]`。
literal anchorはcolumn22で、余分なa倍はない。kappa inverseを用いない。

全149 column DAGと四cached atomを含む162-node計画で

```
T3=B*(kappa dot M15[pair4,original32col3])
det([H27,H31,anchor22,B*col3])=-T3
```

をそのまま評価する。B分母は既存compilerで共通一乗だけ払われている。
全literal分母のQQ/a/B sourceとそのreverse normalizationをcutに保持する。
κ係数は保存された多項式であり、非零や単元性を仮定して逆数を取らない。

元十五行の順は `[0,1,3,6,2,11,12,9,7,13,10,4,5,8,14]`、
十五列は `[22,11,13,25,15,19,14,18,28,30,29,27,31,3,17]`。
既存H2/H3/L/gamma/a row倍率とcolumn3のB一乗だけを使う。
三pair行のanchor消去と多項式Laplaceにより、source domainで

```
A_anchor*D15_3=-T3*D12.
```

ここでD15_3はgenuineなcurrent15全行の選択minorのFULLU既知単元倍である。
GLOBAL full57橋渡しは全六因子環で受理済みで、全Q9 tail/H3/Cdelta効果を保持する。
oldOmega、Delta407/Saux/N6/N7/κ/anchor等のchartは追加しない。
従ってD(D12*T3)はactual full cubic57の十分openである。

## 同じ四次商環で結合積を検査する理由

前と同じe=d=0、prime32003、S=Fp[a]/monic(H_axis)だけを使う。
これはzero-divisorを持ち得る環なので、D12の余りとT3の余りが別々に非零でも十分でない。
必ず結合積 `(D12*T3) mod H_axis` を一回のring積で確認する。
結合積が非零なら、sourceのGauss/LC保持/合法a,B分母制限によりH∤D12*T3となる。
squarefree HとU-coprimeから、同じH成分がD(D12*T3)に会うことが従う。
その二次元openではgenuine D15_3とGLOBAL橋渡しでfull57となる。
H全体を被覆した、あるいは全FULLUを被覆したとは言わない。結合余り零はinconclusive。

literal source cutは9832 column項と4920 cachedκ項を保存するが、値workerは一項ずつstreamする。
全source exponentから有限のLaurent a冪表を作り、aの単元性だけで負冪を処理する。
QQ分母は既存primeで全てunitであることを実行前に確認し、非unitなら候補を変更せず止める。
三次元のT3を展開せず、162個の四次商環ベクトルだけを保持する。
別prime・axis・残りtailへの自動拡張、factor/GCDや旧572再評価を行わない。
fresh ROOT監査と一回GOの前には、この値workerを実行しない。
