# H11 の八 nuisance を含めても二枠の共通欠落曲線は β を優越しない

2026-11-03 / v62。旧 [二枠の被覆](S61_BINARY_HYPERPLANE_FRAME_COVER_20261102.md) は、H8 の五 nuisance の範囲での像次元証明だった。今回、H11 で自由にした三 nuisance も加え、八変数の全値について直接検証した。

S=3e−2 とし、aFC≠0 の共通欠落曲線 B=D0=0 を

    a=−3S(e+1)/50, d=−20/S

と書く。H11 の八 nuisance は

    (u5,u25,u26,u27,u34,u37,u48,u51)

である。上の core 曲線と八 nuisance の九次元 source から W16 への写像を b(e,n₀,…,n₇) と書き、元の primitive 表現の全 Sp4 作用 A(b) を加える。

    M=[A(b) | db/d(e,n₀,…,n₇)]  は16×19行列。

八 nuisance の列は対応する八行の単位行列である。それらを消すと8×11行列になり、50S²を掛けた全係数は QQ[e,n₀,…,n₇] の多項式である。この行列の一行左核を構成し、全11列との積が多項式として恒等零であることを検証した。従って全 source で rank M≤15。

e=2、(n₀,…,n₇)=(1,…,8) では、指定する15次小行列式が QQ で非零である。従って一般微分階数は正確に15であり、標数零の像次元定理から Sp4 で飽和した β 像閉包も15次元である。十分一般の β はこの共通欠落を避ける。

全一次式・分母・左核・非零 minor は [証明書](../../reports/S61_BINARY_NONSQUARE_FRAME_GAP_20261103_v62.json) に保存した。再現は

    .venv/bin/python code/s61_binary_nonsquare_frame_gap_20261103_v62.py \
      --output /tmp/nonsquare-frame-gap-new.json

所要5.14秒。特定 β での微分階数低下を一般化した計算ではない。他の core 曲線、全三次・四次階数、指定点の周囲 t 冪零性を追加主張しない。
