# 五source座標上の有理orbit代表：独立監査

2026-11-04。対象は frozen `Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63`。
既受理の六source、t=2 section、既知固定体の五座標を使う。新しい点・不変式・
source26・stabilizer計算は行わない。

g_u=diag(u,u^-1,1,1)は標準Jを保存する。primitive二形式の元24座標
`x_k dx_i∧dx_j`に対する有限weightを直接用い、全24×16 basis係数で保存した
第一Cartan weightを確認する。元rhoのgenerator1も同じ16×16対角行列となる。
従ってβ_iへの作用は、Lie作用の類推ではなく有限にu^(w1_i)である。

元16有理式へ

    u²=1600r, X=s/(48u), c=xi B u/40,
    v=−32us/h, W=(q+1)/3072, h=s²−25rq²

を代入し、u^(w1_i)倍して保存βbar_iと比較する。作者のquadratic-pair
norm算法を使わず、差の有理式の分母を払い、分子の全monomialを
`u^(2k+e)→(1600r)^k u^e`で独立にfoldする。全16分子のeven/odd部分が零で、
各分母のfoldは非零と確認する。従って二次補助拡大の中の正確な恒等式であり、
保存βbarはuを含まない五変数の有理16-vectorである。非零座標は10。
元t=2のv式も同じfoldで一致する。平方根を有理式と仮定しない。

定義域は表示された非零条件と元source chartの開。特に元の
uv(32u+v)cB≠0を用いる場合、新表示ではh−s≠0も保持する。
βbarの有理式がさらに広い所へ延びても、全境界で同じorbit分類を主張しない。

## sectionのSp4飽和が元像を失わない理由

「sectionと同じorbit」という事実だけから、元六source全体への支配性を
無条件に推論しない。ここでは受理済みのactual raw30 kernel derivation δを使う。
その完全定数体の逆変換によりF6=F5(t)、δはF5を殺す非零f(t)∂tである。
本レビューは旧flow座標のδtをt=2 sectionへ代入し、非零であることを厳密確認する。
これは本稿のshape r,sとは別のflow座標である。

元六sourceのSp4飽和はNEW像Bに稠密である。Bの真のSp4安定閉集合Zの
pullback idealは、actual source kernelと補償群作用によりδで保存される。
元source/Darboux chart上では、その係数はsectionの一般点でregularである。
もしこのsectionがZに含まれれば、generic section点のDVR（uniformizer t−2）
で非零pullback idealは(t−2)^n、n≥1となる。しかしδtがunitなので、
δはその最小vanishing orderを一つ下げ、ideal保存に反する。
従ってsectionはSp4の任意のgeneric安定開を失わず、その飽和はBに稠密。
βbarはそのsectionと二次拡大後の実Sp4変換で同じorbitを持ち、既に有理に
下降しているから、βbarのSp4飽和も同じBに稠密である。

これは五パラメータの有理なorbit代表族である。source-target field index1が
別に証明されれば一般orbit上の有理代表選択として使えるが、このpacketは
その指数やgeneric stabilizer自明性、15次元B自体の純有理性を証明しない。
全25次元Poisson成分の有理性も含めない。

既受理の全次数scale則から、任意のI14について
`I14(βbar)=xi^7 h14(shape)`であり、h14はβbarのxi=1への直接代入から
二次補助根なしで計算できる。これは計算recipeであって、今回I14やh14を
展開・取得したという意味ではない。
