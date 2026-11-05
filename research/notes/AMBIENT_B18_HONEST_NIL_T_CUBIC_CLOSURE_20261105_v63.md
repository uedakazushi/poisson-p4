# 既知三次closureの固定階数とlength completion

2026-11-05。独立査読済み実HF553とsource-only preflightの入力に対し、承認された一回のQQ階数計算を行う。
600×1284の全既知三次sourceは、exact24×1308 presentationへ縮約されている。
固定minor、左右inverse、全1308列factorizationを検証し、factorization streamのhashを保存する。

この階数から得る`24−rank`は、実HF3の上界である。追加の実三次initial relationが無いと仮定しない。
実length16と最初の三層553から残る長さは3。
もし上界が1なら、NakayamaでHF3は零に出来ず1となる。unique-standard-monomialの因子条件で以後も各層は高々1なので、
実Hilbert関数は553111、m6Tnil=0まで確定する。
上界が2または3なら、残るtail候補だけを絞り、正確なHF3を主張しない。

source/gauge/jet係数やGBは新規に計算しない。nil generators・nil cubeは別問題である。

固定計算は1.16秒で完了し、縮約行列の階数は21、既知三次closureの商次元は3であった。
したがって実HF3の上界は3であり、残るtail `3`／`2,1`／`1,1,1`のいずれかをまだ選べない。
この既知closureだけでは三次sourceの追加関係を排除できない。
