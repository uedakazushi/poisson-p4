# 実 B11 上の150係数 gauge 補正と native S2

## 実 source と今回の範囲

必要二次 monomial `r²,z0r,z1r,z4r,z5r` を保持する11次元商 B は、凍結済み
H20 の実商である。前 packet は全400積を照合し、元105根および164個の
source readback を B へ運んだ。本計算はその保存根を使い、元の物理 source
jet または105単位方程式を新しく求めない。旧 r/g 正則商も再除算しない。

本 packet の実出力は次の二つである。

1. actual proper G8 により定義される `S=Res4(S12_by32(G8))` の B 値。
2. native比較列の有効二次代表 `E=(L0+L1)*S`, `L=Res4(N^t)`。

E は各 native128列の個別の actual X2 であるとは主張しない。凍結済みの
全741零-origin結合、特に620零-linear結合に施した後は、既証の
`L2*S0*c_native=0` により真の二次値と一致する。T の長さ・深さの決定は
この出力の別の後続計算であり、本 packet の結果ではない。

## 有限 source 評価

同じ actual CA-F9 `ell=5(s-r)=0` 座標を使い、B の順序は

```
1,r,r²,z0,z1,z4,z5,z0r,z1r,z4r,z5r
```

である。元の二つの48 source の全150 physical tracepair を、同じ scalar
DAG の保存105根から B 上で評価する。全定数・一次は前の native128 packet
の対応成分と一致する。actual monic nuA,nuB も H20→B の同じ写像から得る。

gauge 方程式では補助 etaA を使う rank2 algebra
`B[etaA]/(etaA²-nuA)` を用いる。先に physical tauA と同一視しない。
既存の15 even /15 odd scalar gauge 函数の定数・一次を入れ、全15 constraints
の定数・一次残差が零であると確認する。残る五二次 monomial に対する
30本の scalar residual を保存する。

補正は150 QQ係数であり、既存の QQ origin J30逆を一度掛けるだけで求まる。
J30 は nil中心における非零 mixed block A1 を含む既証の30対30行列で、
determinant は4/9である。補正が `m_B²`、`m_B*m_B²=0` なので未知の
base高次Jacobian項は寄与しない。補正後の **全30 scalar B方程式**が厳密に
零であることを確認する。

## 実 proper frame と G8

同じ scalar gauge を
`C_B=B[tauA,tauB]/(tauA²-nuA,tauB²-nuB)` へ埋める。この algebra は44次元
で、基底は各 B基底に `(1,tauA,tauB,tauA*tauB)` を掛けた順である。

元の source frame は `p=tau+s/2`、実 gamma 群変数と凍結 core 行列から
再構築する。CA frame の正規化が恒等であり、十の固定 fibre 座標が中心値
になることを確認する。同じ physical CA/CB tensor と CB flag から actual14、
retained proper kappa9、および物理F12を得る。その全定数・一次は前の
保存 packet と一致する。

G8 は source-bound proper elimination に従い、五つの単位根からの差と
三つの実 core CI generator の値として直接評価する。反転するのは
`D=1-d²` と `P=5+d²` など既知の単位のみである。全G8の定数・一次も前の
保存値と一致する。G8の原点値は tauA*tauB の非零係数を持ち、零とは
扱わない。N2、R12二次行列、T二次行列は構成しない。

## scalar L0+L1 と出力規約

L0+L1 は、前に保存した actual N の C36値と monic一次係数から **scalar
48対48行列として** Res4 を取って作る。その二次係数を全部除く。native N
係数だけを一次で切り、二次 monicを使って Res4 を取る操作とは異なる。
後者では一部の scalar L2 が残るからである。

actual G8 から28 Koszul +4 free-zero の全source32列を作り、actual monicを
使って Res4 し、48対128の S_B を得る。E=(L0+L1)S_B を計算する。E の全
定数・一次は前の native128 scalar 出力と一致する。scalar L2 の項を
捨ててよい理由は、別の source固有の直和証明で全零-origin結合に対して
既に確認しており、E の各列を真の X2 と誤認しない。

出力 `actual-native128-effective-S2-B11.json` の主要 key は

```
source_columns128_rows48_sparse_B11_effective_second
actual_S_source_columns128_rows48_sparse_B11
```

で、各列は `[row,[[Bindex,QQ],...]]` の sparse配列である。必要二次indexは
`[2,7,8,9,10]`、元の八変数 monomial対は
`[(7,7),(0,7),(1,7),(4,7),(5,7)]` である。元full56の底八行は native列で零。
後続の adjoint は必ず保存された source零-origin結合を先に使う。

## 保存と検証

各 phase は source B値、30対5 RHS・補正、actual14/kappa9/G8、最後に S_B
と有効Eを保存する。全lower-order比較と gauge全方程式を残す。ONE worker
は約1.05秒・52MB、保存18129 QQ slot・約195KBで完了し、60秒・2GB・
150000QQ・32MiB guard内だった。rank探索、T長、depth、348行列は計算しない。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_native_S2_B11_20261104_v63.py
```

既存の report/data/code は不変とし、独立監査は保存150補正を代入するだけで
全方程式・G8・Res4を再検証できる。新しい gauge solve を重複実行する必要はない。
