# 二つの例外二次 core での ρ≠0 支持

2026年11月3日。a=1,c=2、h²+2=0 または1000h²+1127=0 の特殊ファイバーで、
ρ=1 の全135残余式を計算した。元モデルの第136行は t の表示であり、
追加方程式としては使わない。

標数32003では最初の二次式は二つの根3083,28920に分解し、第二は h²+9314 として既約。
従って前者を二つの実際の有限素体で別々に調べ、後者はこの既約多項式から作った
native Galois field GF(32003²) で調べた。分解する二次式の商を体と宣言する操作は使わない。

三ケースとも商の次元0、ベクトル空間長4であり、全法線原点の最大イデアル m による
飽和は I:m^∞=(1) となった。従って被約支持は**標準旗の原点一つだけ**である。
元 t の剰余も零。これらはρ=1の切断についての計算である。

ρ の重みは1なので、各二次 core の同じ標数32003の中央環でρ≠0の全支持は標準ρ軸に尽くされる。
この標準軸上の点では b0 b1≠0 であり、別証明の
[全 h に一様な標準旗混合標数開](S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.md)に入る。
従って周囲除去の議論にこの支持計算を使える。切断の t=0 を直接全周囲へ移したのではない。

ρ=0 の新二次軸と旧二次平面はそれぞれ別の入射・横断証明が必要である。
特にこのρ≠0の支持計算だけで全試験直線や一般五点の尽くしを宣言しない。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_marked_sigma_rho_one_quadratics_20261103_v62.py \
  --rows work/marked-resume/work/flag_boundary_geometry_20261101/marked-vertical-symbolic/global-rows \
  --output work/rho-one-quadratic-replay
```

固定二次 core の特殊化を使う bounded GB であり、過去に停止した generic QQ(h) GB の
同条件での再試行ではない。各M2実行上限は180秒、終了値を個別に回収する。
保存結果は `reports/S61_MARKED_SIGMA_RHO_ONE_QUADRATICS_20261103_v62/`。
初回出力を採用するモードも、凍結全係数から再生成した入力の byte SHA と、
個別に観測した終了値0が一致した場合に限る。
