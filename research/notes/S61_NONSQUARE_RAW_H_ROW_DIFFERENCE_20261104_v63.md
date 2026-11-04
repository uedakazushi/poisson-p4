# 固定row2−row4のtailだけを検査

rawH14×5のrow2とrow4はcolumn0が文字通り同じで、他の非零候補はcolumn5だけ。
この一組の定数差だけを取り、残五座標のうちcolumn5以外が全て零であることを
確認する。tailから除く因子はa,d,F,B,C,JとQQ定数だけで、未知因子のfactorや
逆元は導入しない。

前段の固定四係数の定数span rank3から、column0を消す定数方向はこの差だけ
（scalar倍を除く）である。成功した単元tailなら実際のgenuine16 sourceへ
全係数を戻し、pure部分も同じlazy式として保存する。非定数が残る場合は
正確な既知因子と残多項式だけを保存し、全点rank defectや非全射を推論しない。
追加Gaussian・未知因子caseへのlocalization・一般のminor探索は行わない。
