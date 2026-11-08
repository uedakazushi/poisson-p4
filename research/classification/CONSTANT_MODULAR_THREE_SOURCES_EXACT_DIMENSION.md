# 定数modular三つの完全source：飽和像閉包は各既約二十四次元

2026-11-08、PRIMARY E。producer plane_collision_review。新componentは未主張。
入力は `CONSTANT_MODULAR_DEGREE_TWO_WHOLE_REDUCTION.md` とrootの別作者receipt。
残る完全source I/II/IIIのPGL飽和projective像閉包をZ_I,Z_II,Z_IIIとすると、各Zは既約でdimensionは正確に24。
この下界は選択pointの微分rankに依存しない。新算術／operator／pointは実行していない。

## 1. fixed c の全sourceとf strata

modular vectorをc=-∂tへ固定する。残base(u,z,w)でdivY=1、Yf=f、P=∂t∧Y+Jac(f)。
fの加法定数はYf=fで一意であり、PからYとdfが読めるのでfも復元出来る。
従って(f,Y)表示に余分なpresentation fiberは無い。

以下の各f stratumは既約である。

| source | fの全stratum | dim f | 全Y fiber | 和 |
|---|---|---:|---:|---:|
| I | Hessian rank二、radical linear零、critical value零 | 7 | 10 | 17 |
| II | Hessian rank二、radical linear非零 | 9 | 8 | 17 |
| III | Hessian rank一、radical linear非零 | 7 | 10 | 17 |

rank二symmetric3×3行列のstratumはdimension5：kernel lineのchoice2と、そのquotient上の非退化symmetric2-formのchoice3。
Iではlinear項がimage内の二係数であり、critical value零がconstantを一意に決めるので5+2=7。
IIではlinear項三係数のradical成分非零openと任意constant一係数を加えて5+3+1=9。
rank一symmetric行列はdimension3。IIIではlinear項三係数のimage外openと任意constant一を加えて3+3+1=7。
これらは実際のpolynomial f係数の次元で、basis choice又はnormalizer変数を余分に数えていない。

fixed c を保つbase Aff3は、fにdet(G) f(G^{-1}(x-b))というtwisted作用をする。
三stratumはこの実作用の各一軌道である。
通常のaffine quadratic normal formでq、radical-linear、critical valueを直した後、overall非零係数も残さない。
Iではunused uの非零scaleがdet twistを自由に変える。
IIのf=zw+uはdiagonal G=(a,b,c)でa*zw+(bc)*uへ移るので両係数を独立に直せる。
IIIのf=z²+uはunused wのscaleで両係数の共通det factorを直せる。
全逆は実affine逆である。

従って各f fiberは既証明の完全A10/A8/A10と同型で、全Y coefficientsを保持する。
Aff3×固定f fiberから全sourceへのsurjective写像によりsourceは既約。
fへのprojectionの全fiber dimensionも一定なので、全source dimensionは表の17である。
flatness、generic kernelの壁への外挿又は全higher可解性を仮定しない。
元有限sourceの各値が実際にPoissonという必要十分表示を使っている。

## 2. modular vector choiceとH choice

nonzero constant modular vector cのchoiceはC4\{0}、dimension4。
Aff4はこの全choiceへ作用し、固定c fiberは前節のsourceである。
c自体はPのdivergenceから一意に読めるので、全fixed affine chart sourceは既約dimension21。
座標frameのchoiceを追加parameterとして数えたものではない。

任意のhyperplane H=(h=0)を無限遠に取る。actual homogeneous quinticは非零scalarを除き

```
I:   h³*z*w,
II:  h³*(z*w+h*u),
III: h³*(z²+h*u).
```

quadraticのrankはそれぞれ2/4/3で、Iの二simple factorsはhとindependent。
全例で唯一の三重hyperplaneはhであるため、Hはactual Poisson tensorから復元出来る。
従ってH choiceのdimension4を戻すincidence projectionにpositive-dimensional fiberは無い。

固定hに対するh-Casimir homogeneous quadratic liftは一意でもある。
もしEuler差E∧Wが再びh-CasimirならhW_i-x_iW_h=0。
hはx_iとcoprimeなのでW_h=alpha*h、Wi=alpha*x_i、従ってE∧W=0。
ここでh scaleを独立に数えず、Hはdual projective P4上のchoiceとする。

こうして全homogeneous coneの像は既約dimension25。
非零tensorの同時scaleは自由で全sourceを保つので、projectivizationにより

```
dim Z_I=dim Z_II=dim Z_III=24.
```

各parameter零を削って得たchosen openだけの下界ではない。
三sourceの全像はそれぞれこの同じ既約閉包に含まれる。

## 3. 何が判定され、何が未決か

三つの閉包は相異なる。IIのgeneric quadratic rank四、IIIはrank三、Iはrank二である。
rank<=rは閉条件で、tripleHが唯一なためactualDからresidual quadraticもscalarを除き復元可能。
もし二つの既約24閉包が等しければ、generic最大rankが同じでなければならず矛盾する。
しかし、これだけでは各Zがambient成分とは結論しない。
NEW25等の既知dimension25以上の親への真包含はなお可能である。

全sourceに属する一つの事前固定witnessについて、元350×150 projective Jacobi微分のrank>=101が証明されれば、
Euler24とprojective scale一を除いたtangent上界は24。
上で独立に証明した同じZのlocal dimension下界24と合わせ、witnessはsmooth unique24 component上にある。
その成分は閉既約Zを含むのでZそのものに等しく、whole sourceがこのambient成分に入る。
source点ではないkernel sliceのpointをwitnessとしては使わない。

rank<=100ならこのcriterionは不成立であり、新component又はwhole親を推論しない。
別pointのprofile、minor増加、jet列へ自動継続せず、全NEW25/M11等へのactual confluence／parent幾何へ戻る。
新成分候補を採録する場合は、全priorの同dimension成分とのgenericD／closed-factor監査も別に必要である。
このnoteはその全prior監査又はambient算術を完了したものではない。

rootの純独立読解後にのみ、各whole source一つずつの固定finite witness契約を審査する。
全既知CA/CBを含む初回degree<=2 source reduction、旧失敗、旧reports、global、主稿、Gitは変更していない。
