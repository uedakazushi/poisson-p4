# B18次数零C0表示の独立監査

2026-11-04。作者のsolveを再実行せず、元のactual CA係数回路と保存済みunit系から
七つのCA57 prefix入力を再構成する。四十八係数式と九座標式を各nodeで照合し、
既受理のCA57可逆座標の零制限になっていることを確認する。
五つのpure CA関数について全35coordinate quotientの分子・座標・prefix順序を確認する。
隣接prefixは該当座標零で一致するため、これらの商は全次数で正則である。

保存済みjoint105三次解を、独自の十次元切断多項式演算で元CA回路へ代入する。
CA部分はjoint系の先頭57行でありCB未知数を含まないため、これは同じr,s制限である。
g_Aの三次初期形は非零、四つのBbarは三次まで零と確認する。
既知のactual zero-image membership Bbar=g_A Hと合わせ、全四商の原点値は零と決まる。
三次検査自体から未知のmembershipを推定していない。

actual零weight odd idealは(z七本,g_A)であり、七座標を除いた正則二変数環で
g_Aが非零なので八本のregular sequenceである。従ってC0の全kernelは28 Koszulで尽くす。
作者の保存したideal propagation規則を、各和・負号・scalar乗算について独立に照合し、
元Q0全612列のodd係数をこの八本へ持ち上げた全vectorを再構成する。
旧CBのr/gは元のscalarとして保持し、零locusへの代入・再除算を行わない。

これによりD0の正確な九生成・640関係を受理する。同じtrace pairを表す
C0→H10→D4の全列も再構成し、原点の全関係が零、写像rank七を別のQQ消去で確認する。
D0の最小生成数九は従うが、写像の原点kernel二次元からT0の非零性や有限socleを
推論しない。T0のkernel、元B18の深さ、全nilはこの監査では未計算である。

再現入口は `code/ambient_B18_T0_honest_C0_review_20261104_v63.py`。
新しい空の保存先を`--output`で指定する。全証拠と入力hashは同名REVIEWに保存する。

親回路には古いimplicit系由来の重複表示名・定数nodeがあるため、変数の同一性はnode IDで扱う。
初回checkerの「全node key一意」という仮定はこの親形式に合わず停止した。
数学的入力は変更せず、既存のIDとzero/one sentinelを保つreaderへ修正して検証する。
