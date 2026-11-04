# actual augmented25×32：既知FULLU unitだけの一回Gaussian

2026-11-04。入力は凍結済み `S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63`。
ZW射影25×31とfull H7からのpsi列を合わせたactual source-bound25×32である。
full cubic56×90と混同しない。初段input guardは100k terms/32MiB。

一回だけ60秒・2GB・32MiBのknown-unit Gaussianを実行する。
非零有理数×a,d,F,B,C,Jの整数冪として **厳密に確認できる係数だけ** をpivotにする。
まずa,d-monomial候補を見て、なければ既知六因子を厳密に取り除く。
fill見積り column nonzeros*row numerator terms の小さい候補を選ぶ。
未知minor/factorを反転せず、generic GaussianやGBは行わない。
各完成stepのpivot値・unit factorization・全操作前elimination係数・全25×32行列をcheckpointへ保存する。

intermediate各entryのterm上界150k、output各entry100k、原子的checkpoint32MiBを超えれば停止する。
途中stepを新pivotとして数えず、最後に完全保存されたcheckpointのみを有効状態とする。
source/target indexは凍結inputの固定順。選択順minorのdeterminantは、保存したpivotの積である。

全25unit pivotが得られれば、actualaugmentedmatrixはFULLU全点でrank25。
そのときglobal Cdelta unitimageと既知quad product41によりfullcubicrank57がFULLU全体で従う。
primary計算の成功も独立監査完了前には採用しない。
途中でunit entryが無くなること、制約で停止することは、全点rankの否定やgeneric rank低下を意味しない。
既存のnonzeroDelta34からgenericrank25は別に受理済みであり、今回再評価しない。

本jobはH/Quartic full122 preflightとは別である。全quartic239、genericbeta exhaustion、
既知uniformkernel38・残り42は変えない。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_cubic_zw25_units_20261104_v63.py --output work/marked-square-v63/cubic-zw25-unit-REPLAY
```
