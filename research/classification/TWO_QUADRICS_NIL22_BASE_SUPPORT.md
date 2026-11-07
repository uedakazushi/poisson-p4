# nil22-imageの全非零base括弧：重複H又は旧H_QUARTIC、二quadricなら全て平面

2026-11-07。課題E_reduced_nonnormal、及びGの該当因子型。
対象はh-Casimir二次代表で、affine modular場が

`N=y∂x+w∂z`

となる全支持である。これはJordan(2,2)、定数部c∈im Aの型の全点を
平行移動とGL4で標準化したもの。
結論は、base括弧f={y,w}が非零なら、実Pfaffianは重複一次因子を持つか、
既約四次で旧H_QUARTIC_ESSENTIAL5成分に入る、という全源定理である。
従って**被約H+二つの相異なる既約quadric型は必ずf=0で、本質三変数となる。**
f=0での全quartic支持や、三つの未帰属平面conic源の全成分認識を完了したとはしない。

## 1. 元の全27係数源とbinary discriminant

入力は `CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json` の
zero-based row5（型6）。全90係数に `div P=N`, `[N,P]=0` を課した完全27自由係数源であり、
producerをimportしない独立Fraction再構成は既に受理されている。
ここでPf(Q)=0やquartic既約性を先に課さない。

そのbase括弧は正確に

`f={y,w}=−n₂₄ y²/3+n₁₈ yw+n₁₉ w²`。

元の全Jacobi係数の一つが

`−4(3n₁₈²+4n₁₉n₂₄)/3`

である。従ってfのdiscriminantは零であり、非零fはrank一の平方だけとなる。
rank二を標本で除いたのではなく、元の全係数中の恒等式で一様に除外する。

base p=(y,w)とfiber q=(x,z)へ同じGL2行列Aを作用させると、Nは保たれ、
base括弧は `f'(p')=det(A) f(A⁻¹p')` となる。
非零平方 `f=cℓ²` に対しAの第一行をℓ、第二行を独立な一次式に取り、
第二行のscaleで係数を一へ送れる。従って全非零fは

`f=y²`、`n₁₈=n₁₉=0`, `n₂₄=−3`                     (1)

へ帰着する。零fは別に保持する。全体tensorのscaleや追加スペクトルを仮定しない。

## 2. 全線形式と実際の平方式

(1)の全Jacobiから得る一次係数idealは正確に

`n₀=n₁=n₂=0`, `n₉=−2n₁₄/5`,

`n₁₀=(n₂₆−n₁₆)/2`, `n₁₁=−n₁₇/4`。                   (2)

全一次Jacobi行と(2)のrowspaceが等しいことを、独立Fraction行消去で確認した。
従って全源を失わず、非零定数だけを用いてこれらを消せる。
残る実際のJacobi係数には

`−n₁₇²`, `−(25n₁₂n₁₇+8n₁₄²)/10`

がある。被約支持では順にn₁₇=n₁₄=0が従う。
この二つだけがradicalを取る操作であり、exact nil環の一致は主張しない。

残る全係数からn₃,…,n₈の六変数を、それぞれ非零**定数**係数の方程式で
大域的な多項式として解ける。
各実方程式・係数・多項式解は証明書に全て保存した。
残るJacobi係数idealは正確に

`n₁₂(3n₁₆−n₂₆−2)`, `n₁₂(3n₁₃−n₂₂)`                 (3)

である。有限標本、非線形variableの反転、又は未完了GBは使っていない。

## 3. 全二閉源と全補集合

(3)の全被約支持は、次の二つの**閉じた**源の和である。

* A：n₁₂=0。九自由係数。
* B：n₂₆=3n₁₆−2、n₂₂=3n₁₃。八自由係数。

Bでn₁₂を反転してからAを忘れていない。
両源の全六括弧へ戻すと、全四Jacobiは全parameterで恒等零となる。
`T=wx−yz` と書く。
Aの実Pfaffianは正確に

`−y²{ T+n₁₃w+n₁₅wy+(n₁₆+n₂₆)w²/4+n₂₀/2+n₂₁y/3+n₂₃y²/4 }`。

従ってAは重複hyperplane y=0を持つ。零Pfaffianも含めてこの補集合を保持し、
被約H+既約quartic又は被約H+二既約quadricの枝へ数えない。

Bの実PfaffianFは、全自由係数に対して

`F=−y²T+G(y,w)`, `F|y=0=−n₁₂²`                   (4)

を満たす。これは実際の全Pfaffianについての恒等式である。
n₁₂≠0なら、x,zの一次多項式として係数のgcdは一である。
係数体C(y,w)上で非零一次式は既約なので、Gaussの補題によりFはC[x,y,z,w]で既約。
最高次数四の項があるから、斉次化後にh因子が付け加わることもない。
Bのn₁₂=0補集合はAに既に含まれる。

以上で、f≠0の全点は「重複一次因子」又は「既約quartic」のどちらかであり、
相異なる二既約quadricを持つ点は存在しない。

## 4. 八parameter既約quartic全源を旧23成分へ認識する実witness

旧点 `legacy/v06/nilpotent22.tex` §quarticcomponentの
`P_*=h²J+Q_*` を使う。これは既に全周囲の独立監査が完了している点であり、
全150係数のJacobi微分rank102、実vertical/group/scale rank48を持つ。
従って射影周囲はこの点で滑らかな23次元で、唯一の全成分はH_QUARTIC_ESSENTIAL5である。

old coordinatesをx₁,…,x₄とし、次の実GL4変更だけを用いる。

`(x,y,z,w)=(x₄,x₁,−x₂,−x₃)`。

旧modular場 `x₃∂₂+x₁∂₄` はNへ移り、base括弧はy²となる。
この変更後の旧六括弧は、B源の次の八parameter値と完全に一致する。

`n₁₂=−1`, `n₁₃=0`, `n₁₅=n₁₆=1`,

`n₂₀=n₂₁=n₂₃=0`, `n₂₅=4`。

producerも独立reviewerも、旧Q_*の元の六式から実際のpermutation/sign pushforwardを作り、
変換後の全六式をBの全源と照合した。
この点はn₁₂≠0であり、Bの既約quartic開に属す。

Bの八自由parameter開は既約であり、そのPGL飽和像閉包も既約。
その閉包を含む全成分はこの旧点を通るので、唯一成分H_QUARTIC_ESSENTIAL5である。
よって**Bの全既約quartic源が旧成分へ入る**。
全源が一軌道であることや旧ambient rank102の再計算を追加で要求しない。
ここでは全六式の実source包含と既済の唯一成分定理を正確に接続した。

## 5. f=0は全cotangent/plane源である

元27係数源でf=0ならn₁₈=n₁₉=n₂₄=0。
表示を直接見ると、base括弧は零、四つのbase-fiber括弧はy,wだけに依存し、
fiber-fiber括弧はx,zについて一次である。
従ってPfaffianはh,y,wだけの式となり、本質変数数は高々三である。
相異なる二既約quadricの場合は正確に三となる。
`TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` の全七配置・全六磁気源へ接続できる。

従ってJordan(2,2),c∈im Aの**全被約H+二既約quadric支持**は平面七配置に尽くされる。
CPLUS族とcontact3境界はCPLUS_EXT、二固定配置は二つの新23成分へ入る。
contact4_secant_at_contactも全sourceをβ=0へ移して既知有限閉和へ被覆した。
残二固定配置の全成分帰属はまだ未解決である。
これをB又はAの閉包へ自動的に入れたとはしていない。

## 6. Gへの最弱の次元帰結

平面固定配置のdivisor familyは、essential3-planeの選択Gr(3,5)の六次元と
PGL3の高々八次元から、次元高々14である。
CPLUS parameter配置族だけは、これに一parameterを加えて高々15。
固定したDには唯一の一次因子Hがあるのでh-Casimir liftを固定できる。
全Poisson源のfiberは二anchorのGL2 frame四係数と磁気六係数で被覆され、
projective tensor scale一を除けば次元高々九。
fiberへのbase一次式の加算は追加parameterを必要としない。
実際、`z'=z+ℓ(h,x,y)`, `w'=w+m(h,x,y)` により磁気項は
`R'=R+aℓ+bm−U(m)+V(ℓ)` へ変わり、これは既に全六磁気二次式の中にある。
GL2 fiber変更も上記anchor frameの四係数へ吸収される。
従って全PGL飽和源の次元は、固定配置で高々23、parameter配置族で高々24である。
これは未帰属の三配置にも適用する上界であり、成分性を仮定しない。

未知W優越全成分の次元は25,26,27であるという既済のG次元命題から、
このmodular/因子型の全支持は未知優越成分の一般元にならない。
特殊底Eの残二源帰属を終えずにGのこのgeneric排除が得られる。
他modular型の1+2+2、他の被約非正規因子型、全G又は全分類の完了は主張しない。

## 7. 再現と検証の射程

再利用producerは `code/classification_two_quadrics.py --nil22-base` と
`--nil22-rank1-support`。
後者の全閉源・全Pfaffian・旧点一致は0.623秒で完了した。
独立checker `code/classification_nil22_base_review.py` は標準ライブラリFractionのみを使い、
全元Jacobi係数、discriminant、全一次rowspace、実際の二平方式、全六単位消去、
残るexact ideal、両閉源の全JacobiとPfaffian、旧Q_*元式からのpushforwardを確認した（0.081秒）。
producer/CASをimportしていない。

証明書は `CLASSIFICATION_TWO_QUADRICS_NIL22_RANK1_SUPPORT_20261107.json` と
`CLASSIFICATION_TWO_QUADRICS_NIL22_RANK1_REVIEW_20261107.json`。
各専用jobは `scripts/classification.py gate` 通過後に実行した。
旧H_QUARTICの唯一成分 theoremは既存監査を入力とし、新しい算術確認と区別する。
