# 全四次122×64 finite circuit の独立監査

2026-11-03。`S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63` を、実sourceと可逆row操作のrecipeとして受理する。rank64そのものは検証・主張しない。

## 実sourceと追加chart

元57×91圧縮行列をstreamingで読み、非圧縮hashと164333係数項を照合する。16 genuine cubicの全112 heavy係数を各source行番号・列番号・sparse encodingで独立比較する。最初の41product行の全heavy成分が零であることも確認する。三次row-basisおよび354→319行の等号は既存の独立receiptに固定する。

係数環はaFBC-openで、追加chartは選んだheavy7 minor DeltaHと、168×210 pure quadratic-productsの146minor DeltaQを単元とする。どちらも全coreの既知unitとはしない。full cubic57の点別下界やproduct41の点別下界はこのrecipeの前提ではない。

選んだ七G_Eを列vectorとすると、heavy係数Tは行がG_E、列がH*x_iである。`f=(T^-1)*G_E` のheavy係数はI7となる。この向きを確認し、残る九行からそのheavy成分を引けばpure cubic r_jになる。16行の変換はDeltaHを逆にした可逆変換である。

## 49行の28pivotと21commutator

49行f_i*x_jのheavy係数はunordered pair{i,j}に1を持つ定数行列である。i≤jの28行をpivotとし、残り21行との差を取る。元49列を[pivot28,nonpivot21]に並べればrow変換は `[[I28,0],[*,−I21]]` なのでdetは符号、従って可逆である。全21kernelを個別に確認する。

差は `p_j*x_i−p_i*x_j`。固定規則 `H*x_i*x_j→p_i*x_j (i≤j)` に対し、反対向きの選択との差は保持したcommutatorになる。これはweight4自由加群のnormal formであり、全次数のGroebner basisやconfluenceを主張しない。

28heavy列を除いた後の低重み残余はH*q_iの6行、commutator21行、r_j*x_iの63行で計90行である。元六二次式がpure lightであることを全source monomialから確認する。

## 全33 quartic の保存

元weighted chartのequation weightsとselected34からgenuine quartic33行を再列挙し、319 packetと完全一致を確認する。行77と131の実H²係数b77,b131、および保存u,vを読み、`u*b77+v*b131=(50aBC)^2` を独立に照合する。

二行変換 `[[u,v],[-b131,b77]]` のdetはこの既知unitであり、adjugateによる両側inverseもuniversal polynomial ringで確認する。個々のb77、b131を逆にしない。chiはH²係数unit、psiはH²係数零となる。残り31行から各b_k/U倍のchiを引く操作は可逆な三角row操作で、未知b_kをunitにする必要はない。

先にheavy*light²をnormal formにしてもH²係数は変わらない。従って33行はH²の1pivotと32pure quarticに正確に置き換わる。psi一行と他の31original行の全indexを確認し、33本の実spanを保つ。

四次係数のoracleは、受理済み319 packetと同じactual135式、同じselected34 implicit graph、同じweighted degree4成分である。正の重みによる次数4までの有限再帰で各係数が定義される。未定義の新四次式を導入したものではない。大きな係数をここで展開する必要はない。

## 146 pure方向の除去とcokernel

pure quadratic-productsは168行、pure quarticは210列。DeltaQ chartで146角blockをSchur消去する。既存の独立global upper146は147角minorが全係数で恒等零という事実を与えるため、残る22行のSchur像は局所化係数環で恒等零である。標本点だけで22行を捨てない。

全239列からheavy*light²28、H²1、pure146の175unit pivotを消すと64列が残る。非零の残余行は低重み90行とgenuine32行、計122行。したがって追加chart上で

`rank(full weight4)=175+rank(C122)`

が成り立つ。すべて可逆row/column操作と零行の削除なので、cokernelもC122のcokernelと同型である。

低重みproductsはH²方向を持たず、同じ操作で

`rank(low products)=174+rank(C90)`。

よって90×64の全射は十分条件であり、full122と同値ではない。残した32genuine quarticsが不足を補う可能性を保持する。全33genuine quarticを省いた低重み失敗をfull失敗に読み替えない。

DeltaHとDeltaQの全minorを取る有限chart coverも正しい。ただし対象はheavy rank7かつpure rank146の開だけである。新しいrank64、full四次全射、全三次57、横断長65、ambient除去、任意defect曲線の排除は未証明である。

再現（新しい出力先）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_quartic_122_review_20261103_v63.py \
  --output work/marked-square-v63/quartic122-independent-new
```

GB、新しい体、core点のprobe、parametric minorや四次係数の巨大展開は行わない。
