# prePOWER の h22 チャートにおける純三次九行の support 監査

source-only preflight。凍結済み `PREPOWER_HEAVY_Q_CHART` report `45022e18…` の実 H5 チャートを使い、新たな係数値・小行列式・階数・点・素数は計算しない。Q9 の h22 は既存の十四項係数であり、他の八つの Q9 尾も保持する。

三つの literal H-zero 行は compact [0,2,5]、実135ラベル18,39,68である。この三行については既存5×5随伴行列を掛ける必要がない。実 prePOWER pure 行 P と Q9 を用いて

    h22*P_i - P_i[22]*Q9

の十六列を取る。rho は FULLU 単元なのでこの三行から除ける。これは旧四つの正規化因子も POWER207 係数も導入しない。

## 2+3 ブロックによる六つの連成行

選択H5行列の行順は (1,3,6,12,9)、列順は (0,4,1,3,2) で、ブロック下三角である。先頭二行は列0,4だけを持ち、その既存行列式は det2=-(6/5)B p²。後ろの3×3は既に計算済みの407項 det3であり、既存1752項 det5は det2*det3 である。ここでは再計算しない。

全七 B 行の Cramer 座標は既存報告に5469項として保存されている。四つの B-kernel 行は先に

    det3*G_u - sum L_ui*G_Bi

を作り、残る二つの pair H 成分を2×2随伴行列で消す。二つの pair-only 行 [4,7] は直接2×2消去する。従って六行は新たに5×5随伴行列を展開せずに組める。局所化は既存 det2*det3 と h22 のチャートであり、その補集合を除去しない。

## support 上界と次の実計算候補

全 scalar DAG の入力は既存 checkpoint、accepted_pair_second、新しい大域 Q 完成の小原子だけである。最初の試みでは分母を a の単項式と非零 QQ だけと仮定して source gate で停止した。この失敗は `first-attempt-frontier.json` と元コードに保存する。実際には必要な checkpoint 原子99個に QQ*B 分母がある。B の既存六係数と全係数を直接照合し、以後は a の Laurent 冪と既知単元 B の分母冪を追跡する。一般 gcd や新しい多項式係数は計算しない。

係数を掛けず、Laurent a の指数箱と畳み込み項数の最小値を伝播して support の上界を計算する。和では B 分母の冪ごとの最大値で分母を揃えた上界を取る。打消しは仮定せず、literal zero、同じ根の差、受理済み九つの実 h_i の alias だけを使う。

連成六行の rho 分母は既知 FULLU 因子の指数ごとの最大値で共通化する。これは一般多変数 gcd ではなく、保存された単元 monomial の LCM である。必要な a の冪も明示する。三純行の48係数、六連成行の96係数の上界と、三純行から選ぶ安価な固定3列候補を報告する。候補の structural matching は非零行列式の証明ではない。

次の係数計算には別の資源割当てが必要である。全144係数の上界が150000を超えれば一括展開しない。安価な3×3ブロックでも実際に零になり得るので、まだ Koszul 単元や full57 を主張しない。選択チャート FULLU∩D(det5*h22) の外、p=0、他の h チャートも残す。

再現:

    .venv/bin/python code/s61_nonsquare_prePOWER_h22_kernel9_support_preflight_20261104_v63.py \
      --output reports/S61_NONSQUARE_PREPOWER_H22_KERNEL9_SUPPORT_PREFLIGHT_20261104_v63

上限60秒・2GB・32MiB。出力の support 上界は予測であり、実係数を計算したという意味ではない。
