# Q1 の指定 F-Hasse 計算：ROOT SOURCE 監査

固定4811係数の本文は、独立整数長除法が保存した元197bit共通分母の表示とbyte単位で一致する。
F=500a+90e²+3e−24 の根分子の l1 は117、分母は500である。
Hasse order 0〜16 の全整数partialは input l1×2^16×500^16 未満、従って463bit未満である。
518bit guard、所有量28708、全積回数上界、入力・出力・失敗prefixを含むRAW予算を独立確認した。

固定workerとdriverの全段cursor、最初の非零層の逆向き再積・再生成、一つのaccumulator、
transactional入力l1、safe help、失敗時の保存経路も点検した。監査自体ではFの値を評価しない。
実際の値と独立検証は別のQ1_F_HASSE_ENTRYおよびQ1_F_ROOT_HORNER_REVIEWにある。

証拠：`reports/S61_NONSQUARE_PREPOWER_Q1_F_HASSE_ROOT_SOURCE_REVIEW_20261105_v63/report.json`。
再現コード：`code/s61_nonsquare_Q1_F_hasse_ROOT_source_review_20261105_v63.py`。
