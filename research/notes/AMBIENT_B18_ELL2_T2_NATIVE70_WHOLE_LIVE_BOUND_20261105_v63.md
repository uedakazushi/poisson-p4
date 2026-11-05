# Native70 の値演算前 live 数の明確化

元の c3c2e9e の一般的な予測式は、入力の最後に付け加えた包絡欄を base 数へ戻していなかった。
最終データ15878を単純に二重化して元の包括 scratch を加えると53911となるため、その一般式は使わない。
これは値計算前に検出した count-only 問題であり、元の source code/report/evidence は変更しない。

実際の値 worker に必要な同時所有は別の明示分類で数える。
最終入力の二重所有31756、全 parsed helper は6000、符号付き作業 vector8078、出力と serialization3000、gate/peer/custody900を確保する。
合計49734は元の cap51263を下回る。元の384 bit以上への cap変更も、応答値の再演算もない。
metadata900と最終所有は worker が値演算前・終了時にも検査する。
