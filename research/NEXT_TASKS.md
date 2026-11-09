# 現在の次の作業 — 公表前の独立査読

2026-11-09 / v73-prepublication-review

研究方針を最初に読み、具体的な全検査は research/publication/REFEREE_PLAN.md に従う。

1. CURRENT.md と REFEREE_REPORT.md の今回の結果を読む。A1/A2/A3/B1の独立導出・恒等式を最初から繰り返さない。
2. Library binary downloadが利用可能になったら、Page上の実五volumeを取得する。現状のツールでは実readが拒否され、browser操作ツールも無い。
3. 実五volumeから隔離復元し、SHA・相対参照と四稿の二回組版を行う。保存SHAが一致するローカル巻の検査とは区別する。
4. 査読の保証範囲と全依存を最終確認し、C1が閉じてから公表候補版・新current bundleを確定してmain/Libraryへ保存する。
5. 静的パス索引を完全な動的入力閉鎖と呼ばない。再現で必要な入力不足を発見したら正確なproducer/inputへ戻り修復する。

新しい成分・一般補題・全軌道の探索ではなく、証明上の高リスクを独立に閉じる。
単なるPASS、数式の標本計算、GitHub/Library保存は査読完了ではない。
重大な誤りが実際に見つかった場合だけ、依存する分類定理・台帳・原稿を正しく改訂する。
