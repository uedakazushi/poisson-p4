# V5：六つの関係と全Uの自由57 row module

2026-11-03。V3/V4の保存物は変更していない。今回の公開済み証明書は `S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62`。元63×91matrixのrow moduleが U={aFBC≠0} 上で自由57となること、全点での階数上界57、一般点での階数57を証明した。**全点での階数下界57や自由cokernelは証明していない。** 四次の全被覆・残るbad coreの零次元化も別課題。

公開producerは既存の全135式由来matrixを凍結archiveから取り出してhashを照合し、四単元a,F,B,Cだけの行正規化とその逆の再構成、六つのQQ多項式恒等式を全91列で検査する。正規化にJ,dの因子は不要だった。新しい57minorは (a,e,d)=(1,2,3)、法32003で27632。

五つのgenuine cubic row48,55,59,57,61を先に、最後にquadric-product row14を解く。六対角係数は1,1,1,a,−a,12a/5、三角行列式−12a³/5である。d≠0では最後の14を22へ替えられ、係数−ad/10、三角行列式a³d/10。両方の57×91matrixを公開reportへ保存した。行番号は全て原matrixの0始まり。

## 探索と独立検証

V4で得た三行の削除に加えて、P=a*dの係数座標を使った。全matrixに共通のa⁴を掛けて分母を払い、rowごとの追加正規化をしない表示を用いる。nativeの第三slot名dはこの表示ではPを表す。

core係数の次数をrowごとに限定した有限Macaulay係数行列を作り、2^61−1でsparse echelonを行う。候補kernelをrational reconstructionし、元matrix全列でのQQ恒等式とfree-coordinate independenceを検査する。有限体rank下界とQQ kernel上界のrank sandwichにより、成功した**指定した有限係数空間**のkernelを証明する。

- product係数次数≤7、genuine係数次数≤3：5300未知数、rank5298、二つのQQ関係。原row59=1、57=aなので二行を削除。
- その後、product≤7/genuine≤4は5515未知数でkernel0。product≤8/genuine≤4も7360未知数でkernel0。
- product≤8/genuine≤5：7717未知数、rank7716、一つのQQ関係。原row61の係数−aなので最後の一般依存関係を解消した。

これらのP係数を P=ad として元のQQ[a,e,d]へ戻した。公開検証器は探索時の有限体結果を前提にせず、六つの恒等式を元63行から直接再検査する。source/data/note/reportの現在hashをmanifestに固定した。

`numeric-support.py` の三点のsupportは、使う可能性があった探索用提案であり、零座標が恒等零という証明ではない。実際に採用した最後の関係はその提案を使う前の全row-bound searchで見つかった。公開証明の前提ではない。

## 次の実作業

独立review後、公開reportの `alternative-ad.txt.gz` を出発点とする。この57行表示は、記録した30の単元pivot行を全て保持する。旧unit-spanの八つの追加pivotはcolumnの定数線形結合であり、隠れたrow混合ではない。先行21pivot、追加八pivot、row23の厳密な短いwitnessを対応づけて再現できれば、三次rank欠損のfull-row27×61表示となる。

保持すべき30の原pivot行は

    0,1,2,3,4,5,6,7,9,13,14,15,16,17,18,19,20,21,23,27,28,30,31,35,36,37,38,39,40,41.

row22はこの集合に入らない。旧33行表示に見えた零行は原row48のgenuine cubic重複であり、quadric-product syzygyの零行と混同してはいけない。新しい57行表示に自動的な零行があるとは主張しない。27×61の構成・QQ coker支持・真正四次の全core被覆は、この保存時にはまだ実行していない。

このbundleには、小さいsource/実行ログ/実際のexit0receipt/発見した有限係数を収録する。大きい導出matrixはhash付きで省略し、元入力と変換式・P-chart metadataから再生成できる。証明の軽量replayにはkernel探索をやり直す必要はない。保存時点で自分のsolver processは全て終了しておりlive PIDはない。
