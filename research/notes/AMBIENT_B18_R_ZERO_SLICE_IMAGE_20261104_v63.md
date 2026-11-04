# Actual R曲線の零次数sliceと、同R0内のDVR image

2026-11-04。凍結済みのactual R曲線を、現在のcommon G19切断へ正確に戻す。
raw曲線が既にslice内にあるとは仮定しない。
finite rational group mapと一次unit matrixだけを使い、高次source/group seriesを解かない。

## 1. Raw曲線はtorus-fixedだが切断には入らない

同じactual P0を通る曲線は

```
P(T)=e_y wedge (x²-13T hx+25T²xz, 12T x², 0,
                 hx-13T xz, h²+z²+2T xw).
```

全R membership、Poisson恒等式、Rのこの点でのsmoothnessは既受理packetを使う。
physical variationsのrows16,22,65,112,128は全てactual torus weight0。
全非零constant coefficientもweight0であり、raw曲線全体が固定される。
w=row140はraw曲線で恒等零。

しかしcommon cuts19にはrows16,22,65が含まれるので、raw曲線は現在のsliceには入らない。
sourceparameter Tとgroup19変数gammaに対し

```
(1+gamma0)*G_*P(T), G=I+sum gamma_j E_(GL18_j)
```

というliteral有限mapを構成した。G逆行列は原点Iのfixed Gaussian unit pivotで指定。
全150係数と19 actual cut equationを保存する。
新DAGの独立ADでgroup150×19を既凍結orbit matrixと全係数照合し、T微分も照合する。
cut Jacobian determinantは**-32**。従ってgamma(T)はpointed unit Hensel系で一意である。

## 2. 正規化したslice曲線でも全正weight座標は零

actual fixing torusはGL weights(0,0,-1,0,0)、tensor scalar weight1。
全selected E_ijはconjugation eigenvectorであり、finite group section Gはそのweightに同変。
common cutは各々weight-homogeneousで、baseline P0も固定される。
Tはweight0、group19 weightは既受理の15個の0と4個の1に一致する。

cut解の一意性により、非零weightを持つgamma(T)は全次数で零。
残るgroup factorはtorusと可換するので、正規化後の曲線もtorus-fixedである。
従って**全positive-weight physical座標、特にwは正確に零**。
接線のweightから高次を外挿したのではなく、finite mapの同変性とunit解一意性を使う。

## 3. 実際のR0 imageとparameter

同じnormalization mapからR0の九readback rows

```
27,108,109,112,113,117,118,119,61
```

を取り出した。19×19 cut matrixでgroupの一次補正を厳密に求め、
readback九成分の実際のslice微分を保存した。
その中に非零成分rho_j'(0)があることを確認する。

このrho_jをparameter Uに取り、19 cutへ `rho_j-U=0`を追加した20変数unit系を保存する。
det20=det19*rho_j'(0)は非零。
従ってsource TもgammaもUの一意なalgebraic power seriesとなる。
他の八R0成分をこの同じunit系のfunction psi_i(U)と書けば、

```
I_R,R0=(rho_i-psi_i(rho_j):i!=j)
```

は実際の八生成regular graph idealで、imageはQQ[[U]]。
これはabstractな一変数normal formやprojectionの非分岐性を仮定したものではない。
source RそのものがsmoothでG19を除くと一変数であり、上記unit20系がその同じbranchを指定する。
証明書に全20式・全readback・八graph式をliteralに保持する。

## 4. R支持はT0内で消えない

実際のRとNのcontactは既独立受理の長さ4なので、このR曲線のgeneric pointはNの外。
群作用とslice normalizationは両componentを保つ。
いまこのR曲線は正次数座標が全て零なので、D0の実際のR branchを与える。

D4はactual w-unit tailで、被約部分はNのみ。
D0→D4のimage B0の被約部分は、Nのzero-degree algebraそのものである。
理由はNがdomainでwが非零であり、Nのdegree0部分はw-unit algebraへ単射に入ること、
およびactual D→Nの全射がhomogeneousでdegree0にも全射を与えることによる。
従ってN外のR generic primeではB0が零となり、

```
(T0)_etaR=(D0)_etaR != 0.
```

T0の非零性をmaprankの差だけから推論したものではなく、実際のbranchと零次数mapで確認した。
R-primary部分が閉点でCMであること、T0全体のcurve support、finite socleの有無はまだ主張しない。

## 5. Conductor-ann試験の限界

source57 first coordinate変換から、このR0曲線のCA/CB F9接線と、g_A3/g_B3へ入れた値も保存する。
非零なzero7 firstか非零なg3があれば、conductorの一部はR上恒等零でないことが直接分かる。
全て零でも、conductorへの全次数包含を意味しない。

既知のgeneric R\Nでindividual CA/CBがsmoothという定理を使う場合、
torus-fixed locusもchar0でsmoothなので、対応するzero algebra C0もgeneric Rでsmooth。
よってそのnormalization conductor I0_A/I0_Bはそこでunitとなり、
その冪だけでR支持のT0全体をannihilateすることはできない。
今回、多数のf_i*L1 identityを仮定したり、conductor-ann certificateを主張したりしない。
R-primary部分を分離した後の残りに対するann boundは別問題である。

Code `code/ambient_B18_R_zero_slice_image_20261104_v63.py`。
Report/folder `AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63`。
新しいsource jets、群高次展開、point、prime、GBは無い。
