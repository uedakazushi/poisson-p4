# J=0 full cubic coker: 二つの bounded method と現在の有限問題

## 確定している数学

`S61_NONSQUARE_J0_DEGREE1_UNIT`、`S61_NONSQUARE_J0_DIXON_UNITS`、`S61_NONSQUARE_J0_PRECONDITIONING` の各20261103_v63 packet により、Jを逆転しない実際の cubic cokernel に39回の指定単位 pivot があり、同じ cokernel の18×52行列を保存した。環は

\[
 A=\mathbf Q[a,e,1/(aF B_s C(3e+2))],\quad J=0,
\]

\(d=2(3e+2)/(5a)\)、\(F=500a+90e^2+3e-24\)、\(C=100a+18e^2-9e-2\)、\(B_s=300ae+54e^3+9e^2-24e-4\) である。

最良状態は `work/nonsquare-J0-pid-v63/column-units-v1/primitive.txt`、876,541,698 bytes、SHA256 `039ee56e033a9d576aae8e271483f34e10dd111e87f391eb740251dc90c9489a` である。全単位消去と正規化の入力、保存 witness、known-unit 因子、QQ scalar、再生成方法は直前の packet 群にある。

J=0上で generic full cubic rank57 は、QQの点 (a,e)=(1,2) における非零整数小行列式で証明した。また、固定した20×54 frame で degree≤1 の多項式係数を用いて一行だけに非零像を持たせる direct lift は全20行で不可能である。後者は一般の高次・局所化係数を排除しない。

## Primitive scalar-row GB

上の20×54 frame の実際の第1行について、entryから指定単位の冪と非零QQ contentのみを除いて、54個の primitive polynomial を作った。元のentryの復元等式を全て検査した。QQ[a,e] の scalar ideal のGBと全GB生成元の lift は、60CPU秒/1.8GBの bounded M2 1.24.11実行で完了した。4個のGB生成元に対する全373 lift 項を保存し、QQ sparse多項式演算で **元の行へ戻した四つの完全な等式** を検査した。

うち二つの右辺は `(3e+2)^3` と `a*(3e+2)^2` で、各々14個の非零係数、34項、次数2で表示できる。これは既に消去した行の別の witness であり、新しい pivot 数や full cokernel の結論ではない。係数の分子・分母は最大4294 bitとなり、既存のDixon witnessの1305 bitより悪かったため、新しい大展開には使用しなかった。

最初の Python 後処理は SymPy QQ 係数を `Fraction` に直接渡した型変換で失敗した。M2のGBと全liftはその前に保存済みだったので、文字列を経由する正しい変換に修正し、**保存済み GB を読み直して検算しただけ** で、同じ代数計算を再実行していない。この失敗は元の等式の成立を否定するものではないが、最初の未完了後処理を成功した検算として数えてもいない。

## 二変数 module GB

有理関数体の分母成長を避ける別の方法として、J=0代入直後のcompactな36×70実行列（108,530項、SHA256 `4e8ed3e3c8eb89b614b15ab9a11c91ff51b90e5ef6cca826a27043eb6167e683`）を QQ[a,e] の多項式 module として一度だけGB計算に掛けた。第三の inverse 変数を導入せず、GBが速く終了した場合だけ annihilator と指定単位による saturation を計算する入力である。

この実行は120CPU秒/2GB制限の下で、wall114.5秒後に exit−9 で終了した。ログには入力準備完了しかなく、**module GB、annihilator、saturation のいずれも結果は無い**。終了後の共有cgroup OOM countersは全て0である。終了原因を独立に特定していないので、SIGKILLという観測と実際の制限だけを記録し、未完了を非単位性や正次元支持の証拠にしない。

初回のadapterは、plain matrix のheader `36 70 0` にmetadata二行が無いことを取り違えて、M2起動前に停止した。この読み取りだけを修正して上の一回を実行した。solverを重複実行したものではない。

## 再開する有限問題

未解決なのは、保存した **18×52 cokernelのA上の支持に曲線成分が無い** ことである。generic rank57だけではこれを与えない。十分な具体的証明は、例えば次のいずれかである。

1. QQ(a)[e] と QQ(e)[a] の両方で、指定単位以外は係数関数体の分母として明記したうえで cokernel が零であることを示す。両座標の annihilator があれば曲線支持が排除される。
2. 実際の最大 minor を有限個選び、そのQQ[a,e]上の共通因子が指定単位の積だけであることを exact に示す。
3. 同等の module annihilator/support 証明を作り、元の実行列に対する多項式cofactor等式と、指定単位によるsaturationを明示する。

同じSmith、全逆行列、二変数module GB、大きい展開pivotを無変更で繰り返す再開は勧めない。元の低次数36×70表示、product41の既知の有限欠損証明、genuine16方向の分離など、表示の構造を使う余地が残っている。product41の証明だけでfull cubic57やquartic/ambient nilpotenceが従うとは解釈しない。

入力script、全scalar lift、bounded実行ログ、失敗を含むreceiptを `reports/S61_NONSQUARE_J0_METHOD_FRONTIER_20261103_v63/receipts.tar.gz` に保存した。manifestは各memberのSHA256を持つ。巨大な再生成可能行列は含めていない。現在、子solver processは残っていない。
