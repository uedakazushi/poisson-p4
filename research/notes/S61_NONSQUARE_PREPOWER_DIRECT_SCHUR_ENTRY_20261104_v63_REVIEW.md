# pre-POWER direct P4＋H5 chart と一成分の独立監査

**状態：一成分の算術監査は時間ガードで停止し、未完了。**
下記の `P/a^k` 演算器＋固定 bordered5 determinant は48秒ガードに達した。
各算術段階の完了 receipt を保存する前だったため、7022係数一致や部分的な
数値恒等式も受理しない。上限を拡張した再実行は行わない。
入力・失敗時コード・実行ファイル・ログは ignored work に保存した。
通過済みの型・source・chart 証明だけは、別の
`S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63_REVIEW` に分離する。
以下の数値チェック記述は実装した検査条件であり、PASSした結論ではない。

2026-11-04。対象は `DIRECT_P4_H5_SCHUR` と `U4_T15_SCHUR_ENTRY` の
二つの frozen primary。producer の実行・import、新しい成分・小行列・点・素数・
Gaussian は使わない。新たに評価する対象は指定された同じ768演算と一成分だけ。

## 型付き chart と actual source

全151節点の型、外部ファイル SHA、JSON pointer、全行列積・block・選択・順列・
逆元の形を独立に確認する。P4 は current15 rows `[0,1,3,6]`、original32
targets `[22,11,13,25]`。H5 は compact14 rows `[1,3,6,12,9]`、すなわち
current15 `[2,4,7,13,10]`、heavy targets `[27,31,28,30,29]`。
actual P4 の heavy5 は0であり、全70 heavy 成分の source trace が保たれている。

H5に選ばれた5 source行だけを rho で割り、最後の heavy column31 だけを gamma で
割る。P4 は未変更なので、受理済みの literal determinant
`D4=a^4 C^3 R434` を使える。h22 は逆にしない。H5 の2＋3 block順と既存 minor順も
一致する。`Hnum=[[δ3 Adj2,0],[-Adj3 C32 Adj2,δ2 Adj3]]` の両側積は
`δ2δ3 I5=Delta5 I5` となり、H5 inverse の新 chart逆元は Delta5 だけでよい。

したがって `FULLU∩D(R434 Delta5)` 上、実9 pivot block
`K=[[P4,0],[C54,H5]]` の逆と表示済み L/U により
`L M U=diag(I9,D-C K^-1 B)`。左右の逆も明示され、元15×22と残6×13は同一 cokernel。
未知 projective kernel の自由基底を仮定する議論ではない。
全8個の Q9 の追加 pure tail は実 root のまま残13列に含まれる。
元25 sourceへ戻す covector は、既存の forward Gaussian source operators、
actual S17/S16 と削除行への零挿入の順序付き積を掛ける。transposeや無関係な逆は使わない。

## 指定一成分だけの算術

一成分は compact u4/current15 row5、original32 target15、残行列位置 `(0,9)`。
actual compact H のこの行は pair-only なので、下3×3 cofactor はこの成分に入らない。
独立 checker は全768演算の原子を actual Q2L source、受理済み108値、h9 source、
raw H5、既知 named units へ戻し、各 QQ 分母と a 冪を確認する。
今回の closure に B の分母や h22 の逆はない。

別の FLINT `P/a^k` 演算器を用いる。producer の min-a 指数正規化やコードは使わず、
正の指数 shift と共通 a 分母だけで計算し、last-use後に多項式を解放する。
実 live係数の上限150000を各演算で検査する。実行ファイルと再現入力は ignored work に置き、
巨大な結果多項式を重複して commit しない。

actual pair2 から3 source cofactorsを独立に作り、全5個の Bhat 成分を実 source 値から
再構築する。最終分子については producer の Adj4 contraction と異なり、固定された
bordered5行列 `[[P4,target15],[Bhat_piv4,Bhat_15]]` の120 permutation determinant を使う。
これが保存された Saux の全7022 QQ係数と一致することを必要条件にする。

実 hybrid entry への係数は
`rho_u/(delta2*rhoLCM*D4)`。u行はH5選択5行の外なので rho_u が残る。
この sourceスケーリングと既知単元を含めて、指定された一成分の意味を受理する。

結論は指定chart上の同一cokernelと、この一成分の厳密非零性に限る。
残6行の階数6、global cubic57、quartic239、R434/Delta5 の零点側の除去、
他成分や新しい因数分解は主張しない。

再現：`.venv/bin/python code/s61_nonsquare_prePOWER_direct_schur_entry_review_20261104_v63.py`
