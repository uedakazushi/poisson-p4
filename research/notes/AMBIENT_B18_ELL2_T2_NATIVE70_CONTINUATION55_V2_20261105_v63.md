# Native55 続行の厳密な source-chain reader

最初の501e0a7d worker は c7146 peer が c3c2 digest を直接含むことを要求し、値演算前に止まった。
peer は正しく4bcf gateと606f旧 source peerを pinしていたが、旧 digest はその旧 peerを介した厳密な transitive chainだった。
元のコード、custody0814、空の数学 prefix と BEFORE_ARITHMETIC_FAILURE を保存した。

V2は新しい peer が literal606f reportを pinすることを確認し、その旧 peer が c3c2 を pinすることを確認する。
4bcf、全 source input、worker/helper/external hash の検査は保つ。
変更はこの metadata reader と report識別子だけで、55応答や rank の算術は byte単位で同じである。
同じ128/385、50884/51263、RAW8305349、60秒/2GBであり、最初の15実値を再演算しない。
初回停止は全G8/55/rank算術の前だったため、V2が最初の実55 jobである。
新しい literal source/code/custody peer を経るまで値演算を始めない。
