# nil22 image：全十一源の pointwise 成分帰属の合成監査

2026-11-08、PRIMARY E_multiple_hyperplane。新しい source/operator/point/jet 計算は無い。
既受理の全 source 分類と各全源の actual membership を合成する純証明監査。
**nil22 image の既受理 multiple-hyperplane 入口の全 rank4 被約支持は、既知閉成分の有限和に被覆される。**
新 component 数は増やさない。一般 nil22 cofactor の別の highest-Q/source class へ広げない。
最後の stable_line_conic は root 独立手査読 PASS と全変数 producer PASS を入力とし、
その正式独立 receipt は別ファイルへ追加される。

## 1. 完全入口と全補集合

`NIL22_IMAGE_SOURCE_REVIEW.md` と V2 exact report の既受理定理を用いる。
定数 modular 部が image 内であるため最初の actual affine translation で零にし、
N=x∂q1+y∂q2 とする。入口が与える全 div/Lie/Q wedge N 方程式190×90の kernel19を使う。
ここから

```
P=F∂q1∧∂q2+sum Mij(p)∂qi∧∂pj,
F=F0+(y+div M2)q1+(-x-div M1)q2,
F0 arbitrary quadratic in x,y (six coefficients),
M=-adj(HessD2+HessD3/2)+p tensor V,
D=D0+D2+D3=detM, D0=det HessD2
```

という全源へ戻る。全 Jacobi の必要十分性は割算をせずに
`adj(M)p=-∇D` で検証済みであり、Hessian nondegenerate を最初に仮定しない。
全 GL2 normalizer は tensor と全 F0 を可逆に同時移送する。
cofactor class の source exhaustion を成分帰属又は dimension bound と取り違えない。

binary D3 の四型は三単根／二重＋単根／三重／零。
三単根は Hessian 四点とその exact GL2 識別で二 stable 源、二重根は四源、
三重根は三源、cubic零は Hessian rank2/1 の二源へ尽くされる。
cubic零かつ Hessian rank0 は D=0 なので rank4 に反し、磁気 F0 は Pf=-detM を変えない。
二重根の v/w 零と非零、三重根の w 非零・w0/v非零・w=v0、
Hessian 判別式の全補集合は V2 report で明示処理済み。
平方／立方根 normalizer は非零 branch の actual GL2 選択であり、零 branch は別行に保持する。

## 2. 十一源の終端は全て pointwise membership

| canonical source | 全自由係数数 | 入る既知閉成分／閉和 | whole-source 根拠 |
|---|---:|---|---|
| stable_line_conic | 6 | T17 | 新 true-period coisotropic curve、全150/350/126、全6逆 |
| stable_three_lines | 6 | T03 | 新 true-period coupled edge curve、全150/350/126、全6逆 |
| double_A | 6 | R4_GA9_CUSPIDAL_PLANE_QUARTIC_24 | a=1/4 regular actual endpoint、全6磁気可逆 |
| double_B | 6 | R4_CPLUS_EXT | ε0 regular curve、全free Saito frame＋全6源 theorem |
| double_C | 6 | R4_DEPENDENT_THREE_HYPERPLANES_QUADRIC_24 | 全8源へのactual binding、全6磁気正逆 |
| double_D | 8 | CN31_25 (R4_MODULAR_NIL31_DOUBLE_DOUBLE_H_25) | 全a,u＋全6を持つ regular punctured source→endpoint |
| triple_A | 6 | R4_NIL22_IMAGE_NODAL_CUBIC_FLEX_23 | 全6源の飽和閉包そのものが既受理23成分 |
| triple_B | 6 | R4_NIL22_IMAGE_NODAL_CUBIC_NODE_TANGENT_23 | 別の全6源飽和閉包が既受理23成分 |
| triple_C | 9 | R4_NIL31_IMAGE_THIRD_25 | 全n,c,d＋全6、PGL/Eulerを戻す regular limit |
| quadratic_A | 7 | MPS40有限閉和 | 全λ symplectic actual curve＋全7源の稠密閉性 |
| quadratic_B | 9 | R4_NIL31_IMAGE_HALF_25 | 全n,c,d＋全6を保つ regular limit |

CN31_25 の略称は既 source receipt のもの。現在の catalogue の dimension/metadata を証明の代わりにしない。
triple_A/B は新しく数え直さず、既受理 source image component theorem を用いる。
全十一を相異なる十一 component と呼ばない。

## 3. 各終端のゼロ・critical 補集合

**Stable 二源。** T03 は全六磁気 map det=-3、T17 は det=-2916 の explicit inverse を持つ。
いずれも磁気係数又は invariant を反転しない。
各固定六係数に非空 generic ε source open があり、full polynomial tensor の端点が元源全体。
T17 では k(0)=2,t(0)=81/2,nu(0)=-1/2、T03 では quadric k(0)=-1。
source period/edge の product compatibility は全 ε で保つ。

**double_A/B/C。** A の cuspidal parent a=1/4 境界は source frame denominator の単元内で全6逆を持つ。
B は ε2 を endpoint に使わず ε0 正則、generic ε0/2以外の shared-tangent whole free-frame 源に一致する。
u=1/2 の別 chart も既受理 source theorem が保持しており、磁気 six map det16。
C は generic ε で既知全8源そのものへ入る全 parameter binding で、ε0 は元全6。
三つとも divisor 配置又は dimension/G bound だけの認識ではない。

**double_D, triple_C, quadratic_B。** 全8/9/9源の各点に実正則曲線を与える既独立 receipt を用いる。
parent の単一例が低次元というだけではなく、既全 source の actual punctured PGL 像である。
ε0 が元 source 全点なので endpoint smoothness や local uniqueness を要しない。
自由a,u,n,c,dと全六磁気を一つも捨てない。

**triple_A/B。** 全6 source は irreducible、actual image lower と true-cut support により
その飽和閉包がそれぞれ全 component である既定理を用いる。
全 component の定義が source 全像閉包なので、磁気零点も同じ component 内。
全点が同一 orbit、scheme smooth、component unique とはしない。

**quadratic_A。** αβI!=0 の actual full6 normalizer と全λ symplectic curve を使い、
元 affine A7 が irreducible であることと有限 MPS 閉和の閉性から α=0,β=0,I=0 の全補集合へ戻す。
source射は全域、constant bracket 非零なので射影 tensor は全点非零。
この dense-source closure は元全7源について証明され、chosen magnetic 値の外挿ではない。

全行の D は canonical 非零 D3又はrankS2/1の固定非零項を持ち、全parameterで rank4。
rank≤2の zero-source complement を rank4 点として未処理にしていない。

## 4. 合成結論と射程

任意の既入口の rank4 点を actual affine/GL2 正規化すると、全補集合を含む上表の一行へ入る。
その行は全自由係数について既知 closed PGL-invariant set 内なので、逆の actual 変換も同じ閉和内。
従って入口全支持、全 PGL 飽和及びその閉包が同じ有限閉和に被覆される。

これは E membership の全十一源合成であり、G非優越、source 次元上界又は pointwise divisor label の合成ではない。
一般 nil22 の別 full90 source、shortnil31、Ga9 linear-w、平面1+1+3残源又は全E完了は主張しない。
旧 source review に残る「未帰属六／二源」は当時の歴史的状態として bytes を変更しない。
新しい本稿と各 actual receipt を最終矢印として参照する。
