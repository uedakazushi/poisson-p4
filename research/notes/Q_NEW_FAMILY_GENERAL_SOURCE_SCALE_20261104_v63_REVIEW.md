# 任意次数不変式のsource c重み：独立有限cocharacter監査

2026-11-04。対象は frozen `Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63`。
低七写像の結果からの外挿ではない。保存された全16 primitive有理式と元の
24×16 primitive二形式基底を使い、新しい不変式・点・source26は展開しない。

各非零β座標をQQ(u,v,c,X,W,B)の分子・分母の全monomialとして読み、cの指数が
各々一様で、その差が保存e_iであることを確認する。零座標は零として明示し、
零からeを推測しない。e_iは別に表現のweightとも一致させる。

    e=(-1,-1,0,0,1,0,1,0,1,2,1,2,0,1,0,1).

cocharacter gλ=diag(λ^-1,λ,λ^-1,λ)は元の標準Jを保存する。
元24座標はx_k dx_i∧dx_jなので、その **有限** pullback weightはh_i+h_j+h_k、
h=(-1,1,-1,1)。全24×16 basis係数にこのweightを掛け、保存eに対する
`2e_i−1`が各basis列のweightであることを直接確認する。
さらに標準十Lie生成元のうちindices1,8からhを独立に作り直し、
`−rho[1]−rho[8]=diag(2e−1)`を全16²係数で確認する。単なるweight metadataや
選択した不変式の振る舞いだけを証拠にしない。

従って全source上で

    β(u,v,λ²c,X,W,B)=λ ρ(gλ)β(u,v,c,X,W,B)

が有理恒等式となる。斉次degree DのSp4不変多項式Iへ適用すると、右辺は
λ^D I(β(c))。中心−Iから奇数Dの不変式は零であり、D=2nなら
`I(β(c))=c^n I(β(1))`。sourceではcのLaurent式なので、他のc次数は上の
恒等式により消える。この証明は任意の斉次不変式に適用される。

既受理t=2 sectionのvはX,Wだけの有理式でcに依存しないため、同じ恒等式を
保持する。既知involution固定source体の座標
`K=C(xi,r,s,q,B)`、`xi=40c/(Bu)`ではcの拡大はxiだけを拡大する。
従ってI14=xi^7 h(r,s,q,B)。このhは、元の不変式が既知有限同値を保ち、
I14/xi^7がc-scalingに不変なことからshape体の元である。
将来さらにR,Z0へshapeを変える場合も、cに依存しない有理変換というgateが必要である。

NEW像の既約方程式次数≥26は別の独立受理済み定理である。
従って非零の普遍I14がNEW上で零になることはない。そうなればNEWの高さ1
prime generatorが14次式を割って次数に反する。Inv14が非零であることも
既受理character次元3により分かる。今回その式を新たに展開してはいない。

**条件付きindex1判定。** actual target不変体Lのsourceへの包含L⊂Kは既受理。
全shape体がspectral比とZN等からL内へ回収されたと別に証明できれば、
F4=xi² f(shape)および

    I14/F4³ = xi h(shape)/f(shape)³

からxiもLに入り、L=Kとなる。f,h非零の一般開での関数体の議論である。
ZNのprimitivity・全shape生成・index1という **その前提と結論自体は今回未証明**。
Sp4 ambient不変体の純有理性という別の既証結果も、このsource-targetの指数を
自動的に1にしない。専用receiptはこの条件付き範囲を固定する。
