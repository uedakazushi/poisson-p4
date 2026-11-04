# Honest9 moduleでの初回LOW2負試験：全topが再び二次まで零

2026-11-04。単位Schurで得た正確な9行presentationを使い、最初の低次negative testを行った。
前のambient二行matrixの反復ではなく、全Koszul kernelを保持した9行・2069関係が入力。
actual CA zero7平面へspecializeしてから、二変数(r,s)の二次までの54次元jet spaceで調べた。
specialization後の像に入らないことだけが元のnonmembershipを証明できる。
逆向きの推論はしない。

## 792個の新しいHadamard係数は正確に評価できる

Honest10の新coordinate商は全て元のpure-unit CA source function bのtelescopingである。
actual CA zero7平面上では、先頭七つのprefixは全て原点へ写るので、
七つのz-coordinate係数は正確な定数`partial b/partial z_i(0)`となる。
残り二つは正確に

```
h_r = b(0,...,0,r,0)/r,
h_s = [b(0,...,0,r,s)-b(0,...,0,r,0)]/s.
```

bの原点は0である。従って既保存のsource三次jetからh_r,h_sの二次jetを正確に読み取れる。
この方法で792個全てを評価し、各定数項が既保存の原点値と一致することを確認した。
新しいsource jetやIFT/Newton solveは行わない。

## 古い正則商は未知微分を保持する

古いfull186に含まれるcoordinate商と四つのg商だけは、既知原点値と独立な正order偏差に分けた。
その未知微分を0には置かない。g=rA+sBのA,Bはorder≥2、g tail自身はorder≥3を使う。
単位Schurの10/3 pivot逆元は、これらの偏差を保持したまま二次まで正確な逆元式で展開する。

全2069関係とその全許容(r,s)-多項式倍を含め、不確定な偏差の像を独立に加えて
実際の関係像より大きな有理vector-spaceを作った。
この拡大像にも入らないtopがあれば、実際のT3は非零となる。

## 結果：このspecializationでは検出しない

拡大した54次元space中の関係spanのrankは **40**。
しかし変換後top12の全9成分は、未知偏差を形式変数のまま保持しても、二次まで全て零だった。
従ってnonmembershipを示すcolumnは無い。
T3の零性やactual membership、depth20を意味する結果ではない。
specializationでT3の像が失われ得る点も、そのまま残る。
自動的な次数増加はせず、この試験はここで停止する。

## 保存・実測

Source `code/ambient_B18_T3_honest9_low2_obstruction_20261104_v63.py`。
`reports/AMBIENT_B18_T3_HONEST9_LOW2_OBSTRUCTION_20261104_v63/honest9-low2-certificate.json`に
792個全ての正確な特殊化係数、全2069関係の既知vector・不確定vector、全倍のspan basis、
top12の形式多項式と残差を保存した。
一回3.350159秒、最大RSS375728KiB。60秒/2GB/32MB内。
新しいsource計算、GB、点、p逆元、高次jetは無し。
T3有限長≤36という受理済み結論は維持されるが、その零・非零と全Dのdepthは未決定。
