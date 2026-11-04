# 第三minorのモニック4次ノルム入力preflight

選択済みf=gen10,g=gen15を除き、d次数・項数・総次数・indexの順で
最初のものだけを選ぶ。選択はgen20（88項、d次数4、a次数5、e次数8）。
旧35minorではindex33に対応する。既存Hによるモニック除算でk*d^jの4列を
作り、全係数恒等式を検査する。新ノルムやnested resultantは計算しない。

Hのa分母を払った多項式とkのexact GCDはd次数0なので、提案ノルムKは
QQ(a,e)上で非零になる。しかしKと既存hのcoprimalityは未検査である。
次に小さいノルムを計算しgcd(K,h)=1を証明できれば、a/e射影の候補は有限に
なり、既存モニックHを使ってB-badも有限と分かる。共通因子があった場合は
必要条件の曲線が残るだけで、bad点の存在を意味しない。

有限free112 coverにFULLUの逆元を加えた後の有限性や、special fibreからの
Nakayamaは別条件である。新規点・素数・GB、未知因子の逆数は使わない。
