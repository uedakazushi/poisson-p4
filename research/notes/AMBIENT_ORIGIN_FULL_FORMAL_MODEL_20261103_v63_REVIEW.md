# CA・CB原点の全形式モデルと99関係：独立承認

2026-11-03。対象は `AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63`。
抽象rigidityの実対象への適用と、44+55関係の全kernel性を監査する。
CA・CBを各々別の新座標へ移す主張として承認し、元の二成分の相対埋め込みを同定しない。
全てcomplete local ringの**formal**同型であり、収束する解析同型はここでは主張しない。

## canonical coverのG19 product

受理済みの個別conductorはCM20、三つのsmooth support枝と長(1,1,2)を持つ。
binary認識の全29式の証拠からself枝はCAでU−V/3、CBでU−9V、critical W/Vと独立な接線を持つ。
十九orbit方向に横断するformal group sectionとambient係数空間の補sliceを取る。
action mapのJacobianは可逆で、formal IFTによってG19×sliceがambientへ同型となる。
G不変なactual成分・conductor idealのpullbackはgroup変数によらない。

normalizationとconductor coverはcanonicalなG-equivariant objectである。
有限mapの上でも `(g,e)->g*e` と `(g,e)->g^-1*e` が互いに逆なので、coverも同じproductへ分かれる。
十九parameter sectionがsubgroupである必要はない。R29自体のG安定性も仮定しない。
rank2 finite algebra E/Cのtrace-zero lineはintrinsicなので、そのgenerator変更はunit倍、qの変更はunitの平方倍だけ。
slice上で元のgeneratorを選んでproductへ延ばせば、計算した初期q2をそのまま使える。
他のtrivializationでも非零scalar倍だけであり、rigidityの適用条件は変わらない。

## abstract ring/pair rigidityの全仮定

g19を除いたcoreはCM一次元、edim3、三つのsmooth独立接線枝、generic長(1,1,2)、
全grが `(wr,ws,sr²)`。従って独立のthree-axis rigidity theoremを適用できる。
高次元底上のrank1 torsion-free moduleを無条件に自由とする議論ではなく、product後のDVR証明である。

q=nu mod Iについて、chosen nu-liftの全firstjetがdnuに一致する行列等式を独立に確認する。
従ってq∈m_C²。critical W/V上ではt=0なのでq=0。
actual IFT後のquadratic値はchosen nu-liftの負号であり、
CAではUV/75、CBでは16UV/9。
`s=U,r=U−lambda V` へ変換すると、それぞれ
`(1/25)s(s−r)`、`(16/81)s(s−r)` で非零。
正則化したqは `s² alpha(s)+rs beta(s)`、両定数は非零。
`S=s sqrt(alpha), R=−r beta/sqrt(alpha)` の形式自己同型でqそのものがS(S−R)となる。
勝手なadditive correctionや単なる初期式の置き換えではない。

## regular normalizer内のCohen liftとintrinsic pullback

分類されたsource conductor coverは
`E=C[[g19,w,r,s,t]]/(wr,ws,sr²,t²−s(s−r))`。
q∈m²なので全relationに一次項がなく、edim E=23。
実regular normalizer B29からEへの全射で、Eの23cotangent generatorをBへliftし、
kernel内の六独立linear element z6を補う。Bの29cotangent basisなのでformal IFT/Cohenにより
それらがBの正則parameterとなり、kernelは正確に
`(z6,t²−s(s−r),wr,ws,sr²)` である。

Cのgeneratorはactual Aからliftでき、z6はconductorなのでAに入る。
t²もE像がq∈CなのでAに入る。
従って `nu=t²` として作る
`R=C[[g19,z6,nu,w,r,s]]` はA内の**新しい**正則subring。
B=R directsum R*t、t²=nuである。
任意のnormalization/conductor squareで `A=B times_E C` が成立する：
E像がCにあるbはactual a∈Aとの差がconductorに入り、その差もAに入る。
ゆえに全actual Aは `R directsum I*t`、
`I=(z6,nu−s(s−r),wr,ws,sr²)`。
source imageと任意base changeの交換も、旧CA/CB共通Rとの座標alignmentも使わない。

## 全99 kernel

独立の整数word算術を用い、R⊕R*t、t²=nuのpair乗法で全44syzygyと55pair productを評価して零を確認する。
core Hilbert–Burch matrixは

    [[s,rs],[-r,0],[0,-w]]

であり、符号付きmaximal minorはwr,ws,sr²。
従ってcore syzygyは `sY_wr−rY_ws`、`rsY_wr−wY_sr²` と正確に一致する。
`F=nu−s(s−r)` はnuとの正則coordinate changeなので、七regular sequence z6,Fとcore resolutionを
tensorできる。homological Betti列は
`(1,10,44,112,182,196,140,64,17,2)`、degree2項の44は
21 Koszul＋21 cross＋2 coreであり、全第一syzygyを生成する。

55productで任意のformal式を `R+sum R*Y_i` へ還元する。
二つのYを除くごとにnu*b_i*b_j∈m_R³となるので、この還元は形式位相で収束する。
BのR-basis(1,t)で像が零なら定数部分は零、残りはIのmodule syzygy。
上の44で全て消えるため、99はpullback零だけでなく**全kernel**である。
別の独立M2全elimination receiptもhash照合するが、再実行はしない。

monomial LCMでconductor coreのprimary分解
`(wr,ws,sr²)=(r,s) intersect (w,s) intersect (w,r²)` も確認する。
`U=s,V=r−s` でcritical quotientへqを加えると、両向きの明示係数により
`(WU,WV,UV,U³)` となり、既知source critical模型と整合する。

既受理のcone dim29/edim39/depth21/e2、projective dim28/edim38/depth20/e2を保つ。
この全形式モデルは指定原点の個々のCA・CBについてのもの。
全global分類、元の相対CA/CB配置、actual D、full Poisson schemeへ拡張しない。

独立コードとhash証拠は `code/ambient_origin_full_formal_model_review_20261103_v63.py` と
`reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW/` に保存する。
