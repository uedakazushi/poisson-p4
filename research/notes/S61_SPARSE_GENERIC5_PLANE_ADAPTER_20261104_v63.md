# 疎な一般五変数入力上の既知3平面と42normal座標

元の全134式を残す前段として、普遍トーリックモデルのL51×20、N45×20と二つの楕円切断を、正確な歴史S61五変数入力上で扱う。パラメータはH,L,s0,q0,t0、未知qは45個である。

まず既存の18行・18列selectorの行列式をQQ多項式として求める。零ならそれを保存し、同じ一回のジョブ内で51×20全体の決定的なfraction-free pivotから別の18行・18列を選ぶ。点や素数による選択は行わない。

選択18角行列式をDeltaとし、自由2列を単位行列にした核をT=Tnum/Deltaとして求める。Bareissの各除算を正確な多項式除算と再乗算で検証し、全51×2のL*Tnum=0および自由部分Tnum_free=Delta I2を確認する。Delta非零によるrank>=18と独立な核2本によるrank<=18の両方から、階数18が従う。

一つの既知アフィン3平面は

`q=center+[N*T,Q52-center]*(r,s,sigma)^t`

である。保存されたq行0,1,2の3角式も実際に計算し、零なら適応した3行を選ぶ。多項式列V=[N*Tnum,Q52-center]の選択3行blockをB、その行列式をGamma3とすると、平面の枠はE=V adj(B)/Gamma3に正規化できる。最初の二列のDelta分母はこの比で相殺する。B adj(B)=adj(B) B=Gamma3 I3とE_pivot=I3を検証する。

これにより全45個の未知数の可逆なアフィン座標変換

```
q=center+E*z+std_nonpivot*n42,
z=q_pivot-center_pivot
```

を得る。これは一般係数体でのDelta,Gamma3の局所化だけを使い、未知qの式を単元にしない。全パラメータ上の非零、全域有限性、既知3+1+1による全ファイバーの網羅は主張しない。

元134式の42normal座標への代入、J0行列への代入、Groebner基底は行わない。成果物には入力L/N/center/楕円切断へのハッシュ、18行列式、核の全分子、45×3枠、42個の残余座標を保存する。

再現：

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_sparse_generic5_plane_adapter_20261104_v63.py --output work/NEW-output --work work/NEW-work
```

総60秒・2GB・32MiBで停止し、旧18角、核、平面、normal42adapterの完成したblockを順に保存する。

## 今回の完了部分と停止境界

一回のジョブは11.428秒で、旧18角行列式と20×2核の二つのblockを完成した。その後、未約分の平面3minor/adjugate側の多項式積が、使用したnative共通helperの一成分150k項という停止条件に達したため終了した。reportの全体 `all_passed` はfalseのままとし、平面rank3やnormal42adapterの完成を主張しない。

正の結果は、保存された旧selectorがこの五変数源で実際に非零であること、およびLの正確な階数18である。Deltaは886項、パラメータ次数25。共通分母Deltaを持つTnum20×2は合計34,673項、最大一成分1,252項で、全51×2の核方程式はQQ上で零、自由二行はDelta I2である。旧18行・列・自由列18,19を変更する必要はなかった。これらの行列と全分子は `old18.json` と `kernel20x2.json` に保存した。

ログの最大19,222項は核検証の中間値であり、完成した核分子の項数ではない。Nの各成分は最大3項なので、N*Tnumは一成分150k未満で覆われる。停止はその後の未約分3minor/adjugate側で生じたと分かる。完成した有理核の共通因子を除く方法や、既知定数単元を優先してLを縮約する別手法により、平面の枠を小さくできる可能性がある。ただし、このpacket内で追加のsolverや再実行は行っていない。
