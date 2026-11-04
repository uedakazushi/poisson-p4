# 六変数 NEW25 source の一意な raw-jet kernel

2026-11-03。前稿の六変数 raw jet を使い、**30×26全行列の一般 rank25、一意な kernel とその非零 source 部分**を有理関数体上で証明した。点評価、GB、追加の高次数不変式評価は使わない。

## 行列と scalar の扱い

前 packet の C₀,L は u,v,X,W,B の多項式行列で、元の source は cC₀,cL である。scalar 列には c∂c を取り、全30行を c で割る。c≠0上で可逆な操作であり、行列は QQ[u,v,X,W,B] 上で定義される。

行順序は値の上三角6係数、その後に linear jet の24係数である。source 列は

\[
 \partial_u,\partial_v,\partial_X,\partial_W,\partial_B,c\partial_c
\]

の六つ。全16個の GL₄ 行列 A に対して

\[
 \delta C_0=AC_0+C_0A^{\mathsf T},
\]
\[
 \delta L=AL+LA^{\mathsf T}-\sum_k(Aq)_k\partial_kL
\]

を加える。全四個の perspective 方向 α は

\[
 \delta C_0=0,\qquad
 \delta L=-2(\alpha\cdot q)C_0
              +q(C_0\alpha)^{\mathsf T}-(C_0\alpha)q^{\mathsf T}
\]

である。従って source6＋GL₄16＋perspective4の全26方向を保持する。

## 一回の有界な関数体計算

全30×26多項式行列を保存した。非零 QQ 定数を pivot とする20回の Schur 消去を行い、全 pivot 行・列と逆向きの kernel 復元式を保持する。この段階は多項式環上で可逆である。残る行列は10×6、全397 monomial 項。

一回の native 計算は QQ(u,v,X,W,B) 上で三角消去し、全ての最終零行を各 entry で検査する。残余 rank は**正確に5**。従って元の rank は20+5=25で、kernel は一次元である。ここで使う後半の非定数 pivot は関数体の単元であり、全 source 上の単元と主張しない。

native は逆代入した kernel の共通分母を払い、多項式 gcd を除いた6成分を出力する。その residual10行との積を全て厳密に0と確認する。Python 側で20個の保存済み定数 pivot を逆に辿って26成分を復元し、全QQ分母を払い、**元の30行との全多項式内積を改めて0と確認**した。

結果の六つの source 成分は零でなく、その u 成分は0である。全 source 成分、16個の compensating GL₄ 成分、4個の perspective 成分を完全に保存した。c 列の係数は dlog c であり、元の dc 成分はそれに c を掛けたものとなる。

一回の native run と復元検査を含めて約0.99秒、上限60秒・2GBで完了した。native は field rank に加えて実際の polynomial kernel を検査しており、rank25という既存の次元数だけからベクトルを推測してはいない。

## source quotient への意味と留保

kernel の source 部分は非零なので、純粋な20 group 方向には kernel がない。従って六変数 source から group nuisance を除いた first-jet 方向は一般に正確に五次元で、その source tangent kernel は得られた一つの有理 vector field の張る直線である。

この vector field と compensating group 成分は、raw value と linear jet を一次で保存する完全な証明書を与える。primitive 正規化と Sp₄ invariant はこの group nuisance を消すので、前稿の1025-node invariant 回路の出力は source vector field によって微分されると零になる。

ただし、本稿の rank/kernel 証明だけで、vector field の有限 flow、五つの有理 first integrals、完全な rational quotient を得たとは主張しない。これらは実際の kernel の式を積分して検証すべき次の問題である。source の境界全体で rank が25であること、kernel がどこでも非零であることも主張しない。

## 保存と再現

* `code/q_new_family_source_kernel_20261103_v63.py`
* `code/rational_function_one_kernel_20261103_v63.cpp`
* `reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63/`

`matrix-certificate.json.gz` は全30×26行列、20個の QQ unit pivot、全10×6残差を保存する。`native-input.txt`、`native-kernel.json`、`native.log` は一回の関数体計算の全入力・出力。`kernel-certificate.json.gz` は全26 polynomial kernel と、全30内積の零検査を保持する。使用する前稿の六変数 raw jet とコードの SHA256 を報告に固定した。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_source_kernel_20261103_v63.py \
  --output work/new25-source-kernel-replay
```

新しい β 点や高次数 fitting を追加していない。再現先は未存在 directory を指定する。
