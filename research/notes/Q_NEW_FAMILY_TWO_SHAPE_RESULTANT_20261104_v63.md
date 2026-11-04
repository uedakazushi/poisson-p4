# 一回のexact s-resultant：47,114項のr射影モデル

2026-11-04。受理済み二shape式に対し、**一回だけ**exact FLINT resultantを計算し、非零の整数係数多項式を保存した。消去した変数はsだけであり、結果はrの射影制約である。**NEW25の定義方程式ではない。** rの境界因子、s回収、第三のN₁₂比による消去、actual target体の指数は未解決で、別の計算予算が必要である。

## 入力のlight content検査

`Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63`の三係数表示から

\[
 E_1=A_1+Z_0C_1-Z_1,\qquad E_2=A_2+Z_0C_2-Z_2
\]

の既存分母を払った式を使用した。この節のEᵢは、実際には分母Fn·r(r−1)²J、Fn²·r(r−1)²Jを掛けたQQ多項式を指す。sourceと二shape式の独立reviewの全hashを確認した。

先に**既知分母r,r−1,J,Fnだけ**を調べた。小さいFnのQQ factorizationは一つのsquarefree因子を返し、その再乗算が元に一致した。各入力を、独立target変数の係数三つ（定数、Z₀、Zᵢ）に分け、これら全てに共通する既知分母因子を完全除算で調べた。どちらの式にも除ける共通contentは無かった。入力式全体の大きな因数分解やgeneric gcdは使っていない。

分母openは元の

\[
 \xi r s Z_0(r-1)T N_q J H_nF_n\ne0
\]

を保持する。既存の特性基底への優越は、Z₀,R,Z₁,Z₂が独立であることを与え、局所化した二shape写像は一般的に有限な函数体拡大を与える。従って一般target上の共通curveはこのopenに存在しない。分母上に固定curveがあれば、その既知irreducible分母因子が独立targetの全係数を割る必要があるが、上のcontent検査で取り除けるものが無いことを確認した。これが非零resultantを期待する入力の幾何的理由であり、零になるchart成分を無条件に無視したものではない。

QQ係数を整数primitive係数へ移す行定数は**200、400**。削除因子は両方とも1で、各入力について

\[
 E_{1,\mathbb Z}=200E_{1,\mathbb Q},\qquad
 E_{2,\mathbb Z}=400E_{2,\mathbb Q}
\]

という全係数等号を検証した。整数入力のs次数は3,5、項数は284,1076。

## 正確な行列式の定義と一回のlibrary呼出し

変数順を**(s,r,R,Z₀,Z₁,Z₂)**と固定した。E₁のs降冪係数4個を5つのshiftした行へ、E₂のs降冪係数6個を3つのshiftした行へ置く8×8 Sylvester行列Sを使う。

\[
 \mathcal R(r,R,Z_0,Z_1,Z_2)
 =\det S=\operatorname{Res}_s(E_{1,\mathbb Z},E_{2,\mathbb Z}).
\]

各matrix entryは、保存した二入力のs係数tableへの参照として完全に定義した。巨大なcofactor bundleは生成していない。QQ入力のresultantとの差は非零定数200⁵400³だけである。

既存CloudのFLINT **3.1.0** の`fmpz_mpoly_resultant`を、整数6変数ctx・消去var=0で**一回だけ**呼び出した。APIのint return値は多項式ではなく成功statusであり、実際にはstatus=1、返った多項式は非零であった。status=0と零多項式を別々に扱う実装とした。別algorithm・別素数・点・source26・r消去への再試行は行っていない。library内部のexact modular algorithmは、外部のβ点評価ではない。

## 保存した展開と範囲

返った整数多項式は

\[
 \#\mathrm{terms}=47114,\qquad
 \deg_{(s,r,R,Z_0,Z_1,Z_2)}=(0,72,28,6,5,3),
\]

最大係数201bit。nativeは約39.27秒、content前処理・compile・保存を含むshared実行は約41.32秒で、120秒/2GB/全artifact32MiB以内に完了した。

canonical係数表は、各項を

```text
decimal_integer_coefficient e_s e_r e_R e_Z0 e_Z1 e_Z2\n
```

として並べた3,008,878 bytesで、SHA256は

```text
9a6b7a4223125bf7f195c8616b963dcd41fb0afe1e8d8f3ee7585504e3d01871
```

である。`resultant.txt`の最初の`6 47114`行を除いたhashで、wrapper側でも全行を再読して一致を確認した。header込みfile hashはreportの証拠一覧に別に保持する。

metadataとして全項の最小指数を読むと(0,22,0,1,0,0)で、**可視のmonomial因子r²²Z₀**を持つ。これを除いた別の展開や、r−1など他の境界因子の因数分解は実行していない。raw次数72をactual一般fiber次数と扱わない。元のopenではrとZ₀は単元なので、この可視因子が零集合を変えないという局所化の注意だけを残す。

この大きなlibrary出力について、別algorithmによる全行列式恒等式監査はまだ行っていない。正確な8×8 determinantの入力・定義と、得た展開を分離して保存し、statusを`nonzero_resultant_saved_pending_independent_identity_audit`とする。raw resultantだけではsの各候補の回収、分母条件、N比との整合、NEW方程式・既約性・target indexを決めない。

## 再現と引継ぎ

producerは`code/q_new_family_two_shape_resultant_20261104_v63.py`と同名C++。

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_two_shape_resultant_20261104_v63.py \
  --output work/two-shape-resultant-replay-report \
  --work work/two-shape-resultant-replay-native
```

新しいoutput/workを指定する。shared120秒にlight content検査とcompileを含み、nativeには残り時間から保存用3秒を引いた時間だけを与える。OpenSSLは既存libraryをSHA256のstream計算だけに使用し、networkは使わない。返った係数表の保存量が残り32MiBを超えればhashと次数だけを残して停止する仕組みだが、今回は全係数が保存できた。

一次保存は`reports/Q_NEW_FAMILY_TWO_SHAPE_RESULTANT_20261104_v63/`。
`content-preflight.json`、`integer-input.txt.gz`、`determinant-recipe.json`、`resultant.txt`、library status・次数・hash、実行frontier、入力・証拠hashを含む。次段階は、独立のidentity audit、明示した局所化を保つ境界因子処理、subresultantによるs回収またはN比とのr消去の、別々に予算を定めた検査である。これらは自動開始していない。
