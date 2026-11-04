# rank4成分と実S61の有理積：独立レビュー

2026-11-04。原著 `S61_RATIONAL_COMPONENT_TRANSFER_20261104_v63` を監査した。
新しい再構成・source jet・評価点は使わず、受理済みの全S61 ideal証明を再利用する。

固定点は最初から `[1:0:0:0:0]`。その四次元接空間のskew値Cの六係数を記録する。
cはこの四次元添字でC03であり、projective五座標のspatial pairでは(1,4)に対応する。
独立にRの行をskew Gram-Schmidtでpairingし、逆は連立式から直接

```text
Rinv=[[1,0,0,0],[e/c,a,1,0],[f/c,b,0,p/c],[0,c,0,0]],
p=af−be+cd
```

と求めた。R*Rinv=I、Rinv*J*Rinvᵀ=C、detR=1/pを全六変数で確認した。
必要な単元はc,pだけで、平方根・有限被覆は不要。

## 実際の旧perspective gauge

affine perspective `phi_alpha(y)=y/(1+alpha*y)` の逆での全微分を独立に計算し、
`h(I−y alphaᵀ)`、h=1−alpha*yを確認した。任意のskew Mに対する
`(I−y alphaᵀ)M(I−alpha yᵀ)=M+E∧(M alpha)` と、radial cubicに対する
余分なhの相殺も全成分で検査した。従ってdegree0..3の有限pushforward式が使える。

線形部分は `L'=L−2(alpha*y)J+E∧Jalpha`。
旧S61の四つのpivotはraw61..64で、実際のresponseは
`diag(−3,−3,−2,−2)`、det36。
従って旧normalizing alphaは `L12` の四係数をそれぞれ3,3,2,2で割ったもの。
primitive trace gaugeの `grad(L01+L23)/5` で置き換えない。
有限変換の逆を使うので、Qの二次shiftも全て保持される。

## 全schemeの積とcomponentへの制限

値Cを記録し、diag(1,R(C))で値をJへ移し、旧gaugeのalphaを記録してNへ正規化する。
逆はNに逆perspective、次いでdiag(1,R(C)^−1)を作用させる。
perspectiveは固定点と値を保ち、記録alphaは定数可逆responseから一意に復元する。
accepted55消去は非零QQ定数のみを割り、154から134への変更は20のideal恒等式。
従ってこれは根基や接空間だけでなく、元Jacobi schemeに対する
`U6×A4×S61` の可逆表示である。

connected PGL5およびtensor Gmは各被約既約componentを保つ。
表示した群変換と逆もそのcomponentを保つため、積は個別componentへ制限できる。
その積が既約ならslice Yも既約である。componentを被約schemeとして扱う場合も、
群作用のfactorizationとfaithfully-flatな積からreduced sliceの主張が正当化される。
無関係な複数S61成分を同一視していない。

一般rank4 componentには、非退化評価点を固定点へ移すPGL変換が存在する。
さらに固定点を保つGL4でcを非零にできるため、このchartはcomponent上で稠密。
評価点の自由度4を次元へ追加する操作ではない。

projectivizationは**正規化前にc=1というtensor-scale chartを取る**。
従ってprojective componentは `U5×A4×Y` と双有理で、次元差は正確に9。
S61上のscalar作用をさらに商にするという推論ではない。
全式はQQ上なので、幾何的componentがQQで定義され、YがQQ上で有理ならXもQQ有理。
C上も同じである。これは `Y_n rational => X_(n+9) rational` という条件付きtransfer。
新成分Y16の有理性、全Poisson25 componentの有理性は本packetでは証明も仮定もしない。
