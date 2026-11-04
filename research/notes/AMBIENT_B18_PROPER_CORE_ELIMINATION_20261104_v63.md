# 十二 flag 商から四変数七多項式へ: 五つの逐次 unit-linear 消去

2026-11-04。元の H0,...,H11 は実 E4=QQ[[source9]]/(d*H12) の凍結済み十二多項式。
Hselected の二つを同時に線形と確認し、残り三つも一変数の unit-linear 式の場合だけ消去する。
GB・局所標準基底・CM・D0 の主張は行わない。

## 1. 二つの線形式

H0,H4 の(v3,v4) coefficient matrix は

    [[0,2(d²+5)], [8(1-d²),4d(-3a+d+v2)]]

で determinant=-16(1-d²)(d²+5)、中心値-80。
従って v4、次に v3 を、この二つの明示単位のみで有理的に解ける。
source九変数と scalar/cone の範囲は変えない。

## 2. 三つの有界検査と逐次 graph ideal

次は H1 を v5、H2 を v7、H3 を v8 について検査する。
各式がその変数に正確に一次で、係数が原点で非零、かつ constant 部分が未消去の
後続 pivot を含まないことを確認する。残り三式に関する次数・単位検査は三回だけである。
どれかが失敗すれば algebraic root を選ばず、小さい暗黙模型を残して停止する。

全て成功した場合、各関係は先行 root を代入した後に

    unit_i*(pivot_i-root_i)=0

となるため、selected五式の ideal は exact rational graph ideal に等しい。
分母の単位値と全五 source 恒等式を保存する。
三次/四次 jet を使って近似 root を作る処理ではない。

## 3. 七つの実残差

最終 base は(a,d,v2,v6)。残る七つ全ての元 Hrow に実五 root を代入する。
分母は D=1-d²,P=d²+5 の明示冪と既知 QQ 定数に保ち、
既知 D/P 因子だけを正確に消去する。未知 polynomial を単位反転しない。
primitive な七整数多項式 Ki と、元 Hrow へ戻す全 unit multiplier を保存する。
同じ formal flag d は free base 変数のままである。

従って core 商は QQ[[a,d,v2,v6]]/(K7) となるが、E4全体は

    QQ[[u1,...,u5,a,d,v2,v6]] / d*(u1,...,u5,K1,...,K7)

である。五つの d*ui 関係を捨てない。
この exact source polynomial input から先の ideal membership や標準基底は、別途許可された後の作業。
この packet は全 core のCM、socle、D0深さや actual-normalizer比較 kernel を決めない。

各積の前に sparse-degree/term bound を確認し、source・既存root・既存残差・現在の積を合わせた
expanded term guard を150000以下に保つ。出力32MiB、実行60秒・2GB。
実際の term 数と unit 値は添付 report に記録する。

再現:

    .venv/bin/python code/ambient_B18_proper_core_elimination_20261104_v63.py --output work/B18-proper-core-replay
