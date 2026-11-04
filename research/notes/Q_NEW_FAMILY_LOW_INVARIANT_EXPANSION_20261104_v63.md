# 五変数 source の低次数 invariant を点なしで完全展開

2026-11-04。事前上界に従い、一回の native FLINT 整数計算で F₄,p₀,…,p₄,N₁₂ の全係数を展開した。**新しい β 点も、24次の点評価も使わない。** 41個の24次 monomial の関係探索は本 packet には含まない。

symbolic t=2 section で c=1 と置く。

\[
 v=-98304X/D,\quad D=147456X^2-(3072W-1)^2.
\]

既に検証済みの β₁₆ に対し b=1080BDβ とすると、全16入力は小さい整数多項式になる。二共変量の有理係数も4で払い、C̃=4C(b)、D̃=4D(b) を作る。

* F₄ は tr(−C̃+2D̃)²、全体の規格化は16(1080BD)⁴。
* pᵢ は det(tC̃+zD̃) の全24置換から同時に作り、規格化は256(1080BD)⁸。
* N₁₂ は整数 word

\[
 8\operatorname{tr}(C\!D^2C\!D^2)
 -4\operatorname{tr}D^2\operatorname{tr}(C^2D^2)
 +(\operatorname{tr}D^2)^2\operatorname{tr}C^2
 -2\operatorname{tr}D^4\operatorname{tr}C^2
\]

から作り、規格化は32768(1080BD)¹²。ここで word 内の C,D は整数化した二共変量である。

tr C²=2tr D²、tr(CD)=−tr D² と F₄²=400p₀+1600p₄ を source の整数多項式として直接検査した。native code は全積・和の項数を検査し、上限を超えれば停止する。

有限 involution と c-weight により全出力は u+X の指数和が偶数である。各 nonzero term でもこの parity を確認し、

\[
 u^iX^j\mapsto a^iY^{(i+j)/2},\qquad a=u/X,\ Y=X^2
\]

で四つの shape 変数 a,Y,W,B へ移す。整数 content と既知の B 単項式因子だけを取り消した。D 因子は取り消さず、分母の D 指数は4,8,12をそのまま保持する。

|出力|項数|係数最大bit|分母のB指数|分母のD指数|
|---|---:|---:|---:|---:|
|F₄|96|127|2|4|
|p₀|336|221|3|8|
|p₁|1040|245|4|8|
|p₂|1061|256|4|8|
|p₃|973|256|4|8|
|p₄|893|256|4|8|
|N₁₂|1700|341|5|12|

各圧縮JSONは numerator の全整数係数と単項式指数、正の整数分母、B,Dの分母指数を持つ。scalar c を戻すと F₄はc²、pᵢはc⁴、N₁₂はc⁶を掛ける。

native の最大 intermediate 項数は4083、全925 polynomial multiplication、実行0.099秒。compile・整合確認を含む wrapper は約1.43秒。全七出力は非圧縮398,534 bytes、圧縮計約176KBで、60秒・2GB・32MiB・output100k/intermediate150k項の全上限内で完了した。native run は一回で、自動 retry は行わない。

これは全低次数 invariant map の完全な係数記録であり、NEW25 の式そのものではない。点の有限 rank を恒等式とみなす操作はない。

保存：`code/q_new_family_low_invariant_expand_20261104_v63.py/.cpp`、
`reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63/` の七つの `.json.gz`、native input、word log、report。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_low_invariant_expand_20261104_v63.py \
  --output work/low-invariant-expansion-replay-report \
  --work work/low-invariant-expansion-replay-work
```

入力には実族から検証済みの β、symbolic section、全二次共変量係数と事前検査の hash を固定する。未存在 directory を指定する。
