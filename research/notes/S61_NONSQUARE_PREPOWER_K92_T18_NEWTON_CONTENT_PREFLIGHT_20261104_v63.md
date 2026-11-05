# 同じK92/18の重み付きNewton・QQ content見積もり

2026-11-04。前の238909格子・64.5MB見積もりは正しい緩い上界として変更しない。
同じ十H3項とa lambda419を全て保ち、parameter係数の積・新値・点・素数を計算しない。

各元source atomについて、九つの非負整数重み（2a+eを含む）の最小最大を元指数から読む。
積は各支持函数の区間和、和は二つの区間の包絡として伝播する。a-Laurent shiftと
既知B分母のclearも同じ規則で追跡する。最終格子数はa,eの二重loopで全重み条件の
d区間を交叉し、点のリストは一切作らない。

同時に、各atomの正QQ contentを元有理数係数から求める。以降は最大contentと主張せず、
既知の下界因子だけを伝播する。積では因子を掛け、和ではQQ gcdを取り、
残る整数係数のl1上界を整数倍率で加える。これはGaussの整性を使うscalar処理であり、
多項式gcdや因数分解ではない。既知のQQ因子とa monomialはFULLU単元なのでsource式を失わない。

元JSONは必要fieldだけを保持し、parse一時peak・元parameter項数・15000作業枠を数える。
同じsource係数列のコピーをcache keyへ蓄積せず、元referenceとQQ単位・a shiftで同定する。
人工小boxだけで区間counterを照合し、実source点の数値検査は行わない。

この報告はsource/metadata gateと出力の新しい上界である。scalar workerの全phase寿命・
scratch・RAWを満たしたという主張とは分ける。失敗または未確定ならそのfrontierを保存し、
大域停止や全FULLU57へ読み替えない。
