# 全 HB 六座標と三 pivot ell 三次：独立 SOURCE・実装監査

独立報告 `89de14dc...` は主 SOURCE `6cbc0f67...`、未実行 worker `54ee6475...`、
custody `4b751cee...` を新しい値を読む前に確認した。
三つの mu03 は row18、60、67 の delta_ell³ であり、既存 lambda_r ell² とは別である。

実 physical の first7 と三つの二次 pivot から H と T を再構成し、
(-B2)T=I3 を確認した。P3_A=T delta3、P3_B=HᵗT(delta3-delta2_mixed beta2) を
全七 source から作り、free四行と graph の既証零を使ってから三 pivot に縮めた。
未知 mu は形式葉のまま、全6×2 affine operator を独立の二変数多項式で再構成した。
負の beta carry、r³ carry、全 r ell³・ell⁴ 座標を保持し、旧行零・三への射影も
数値代入ではなく係数ごとに一致させた。

P0、P1、P2 の実二次は r²H、r ell H、ell²H そのものである。
従って二つの HB の全三次相殺は ell³ 商を使わない。
P4 は r/ell を掛けると五次以上、z 四次の firstframe 補正と H 四次の二次pivot補正も
A/B carry が五次以上になるため、追加 delta4・mixed source3・beta3 は不要である。

7095個の originalID 付き演算の最小 closure、同じ105 unknown/equation、
315個の低次入力と18個の低次 root を再確認した。
正の majorant と独立の閉じた逆数上界により、root 演算前71bit、
affine/rank/両側逆行列検算まで783bitを確認した。
last-use162、値 worker の全所有量は peer5000枠込み58222で、静的75000以下である。
phase1 の root 出力を post 計算前に固定し、post 停止時には root を再実行しない。

私の最初の T 検査は ell² 商 raw33 に存在しない ell² monomial を読もうとして停止した。
そのコードを ATTEMPT1 に保存し、別途検証済みの full physical frame の ell² 係数へ
正しく結び直した。主証拠には問題がなく、変更していない。
主側 custody の旧 path 誤記も算術前の失敗として別に保存されている。

本 PASS で親の条件付き一回実行の前提を満たした。
二 HB が rank二なら全 HF4=4、b=2、e は三又は四である。
rank一なら既知の二 source だけの結果であり、全 HF4 は未決の四又は五のままである。
新しい値はこの SOURCE 監査では評価していない。

再現コードは `ambient_B18_C_full_HB_pivot03_formal_review_20261105_v63.py` と
`ambient_B18_C_full_HB_pivot03_source_finish_review_20261105_v63.py`。
