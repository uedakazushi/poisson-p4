# 十二角軸sourceと固定Bareissの一回限定worker

2026-11-05。源214のsigned axis評価と十二角determinantを一回の予算で行う準備。
旧579整数係数は一切再計算しない。追加23成分だけを元literal sourceから計算する。
入力は67個のliteral atomをstreamでコピーしたもの、214個の有限DAG、旧579、
既存normalization記録である。全checkpoint/巨大q-DAG objectをruntimeで展開しない。

source profileのB(a,0,0)=-40とa単元だけを使い、e,d非負のgeneric入力であることを
全literal atomで確認する。H3のrho13/Lは8000000/16119。追加列にも行一＋二行二を施す。
旧row/column a-dualを保持し、新row/column shiftは0、B clearも0。
row LCDの整数追加倍はmetadata通りにし、全144成分の整数化とsource再読込を検証する。

Fraction source temporary238bit、固定Bareiss temporary1984bitのcompleted-source上界。
実hard guard2048bit、全source準備19212、native33275、RAW閉包22.24MB、45秒soft/60秒hard。
binary M12BP001は全quotient witness、M12AS001は全matrix/buffer/cursor/一時整数を保存する。
fixtureは負数・零・2048bit整数のserializationだけを検査し、source値やpointを使わない。

最終axis determinantが零ならaxis rank十一という結果のみ。generic determinant零ではない。
非零ならgeneric D12 numeratorの非零性と、新旧Omega chartのfull57十分条件に限る。
独立peer前に新しいpivot/global57/fullcurve239/尽くしを主張しない。
