# HEAD29 の coherent 接続式と既存 cache：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`295b19f1a80e7fa512c11218501f32123d8987db2b18b9bf79b84cd527c0b4b7`。
HEAD29 SOURCE `f18d6c07…` と cache availability `fc451d7f…` を固定した。
head29 の応答値や新しい root の係数は評価していない。

旧自由56行の row29 は実 physical63*tauA で、sigma0=e29である。
元 Q48 frameの B 成分全体は厳密に零なので、B sectorに属する
sigma1_z の surviving左 adjoint と L*m2*sigma0 は零になる。
一方、`−L*m1_r*sigma1_r` は残る。後者を省略してはいけない。

独立 reviewerは式の項別転記ではなく、q7 の必要な入力係数から
`V4*sigma0+V3*sigma1+V2*sigma2` を ordinary convolutionで作り直した。
全24個の unprojected sigma1_z/sigma2 source成分を形式変数として残し、
最終 contractionで不要成分が消えることを確認した。
mixed qHzi、B二次 mixed係数、genuine kappa H二次補正も含め、
genuine二座標とcanonical二座標の全四未評価多項式が一致した。

sigma1_r の A curve quotient は delta63_r4、H curve quotientは
`−bB63_r3`、A/Hの座標Hadamard成分は delta63_zi*r と `−bB63_zi*r` である。
H curveは shifted B curveと混同せず、元 `(z7,r²)` sectorを使う。
この headでは新しい mixed third coefficientを必要としない。
低二座標は sourceだけで零である。

raw33から plane9への写像を、全9座標の全33基底への作用297個で独立検証した。
plane9は `QQ[r,ell]/(ell²,(r,ell)^5)`、写像は z7を零へ送る。
P(1)=1と全作用の intertwiningは、cyclic quotient全体で
P(f)=f(0,r,ell)を保証する。単なるmonomial表示の見かけに依存せず、
非斉次 normal-form carryも全て含めた検証である。

したがって実 C2の physical63 delta slot31/32は、literal Taylor係数
`delta63_r4=−1236125/1728`, `delta63_r3ell=−91325/288` である。
これらを「既存入力」という別形式 labelへ付け替えた全四式も一致した。
ここで残る未評価 labelは delta63_r²ell²と、既存 joint105からの
bB63_ell³の compositionだけである。前者の四出力への感度は
`(1,0,5/3,−121/48)`。既存 B6 jobはell³を捨てているため後者を決めない。

reviewerの初回 accessorは全native入力をpureplaneと誤限定し、
必要なmixed H二次項を読む前に停止した。これはauthorが保持していた項であり、
独立 accessorを修正した。次の停止はsingleton callback labelをtupleでなく
stringとして期待した局所 assertionだった。両 code snapshotと失敗記録を保存した。
authorの凍結ファイル・数式は変更していない。

独立SOURCE所有上界80,090、実行0.0563秒、2 GB・55秒内で終了した。
高次functionalの実行コード・正のheight/所有量・custody・新しいROOT GOを
このSOURCEだけから認めるものではない。元 H0・Kbar指数・nil三乗の結論も加えない。
再現コードは `code/ambient_B18_C_connecting_head29_source_review_20261105_v63.py`。
