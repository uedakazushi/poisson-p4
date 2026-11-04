# 独立な六 invariant の24次41部分空間：点なし係数検査は単射

2026-11-04。一回だけ許可された点なし係数検査を行った。F₄,p₁,p₂,p₃,p₄,N₁₂ の41個の24次 monomial の source 制限は **QQ 上独立**である。従ってこの41次元部分空間には NEW25 方程式がない。**全24次不変式空間を除外したのではなく、一般の次数下界は24のまま**である。

## 部分空間と exact source 係数

六 generator の β次数は (4,8,8,8,8,12)。universal rank6 により六つは全W上で代数的独立である。対象は次の全41 monomial である。

* N₁₂を含まない35個：F₄⁶、F₄⁴pᵢ、F₄²pᵢpⱼ、pᵢpⱼpₖ。
* N₁₂を一つ含む5個：N₁₂F₄³、N₁₂F₄pᵢ。
* N₁₂² の一つ。

展開済みの六つの exact source 有理式は、共通分母 B¹²D²⁴で41個の整数多項式にできる。各 generator の正の整数分母の積も列ごとに保存し、最後に universal monomial へ戻す定数 rescaling を失わない。

全 monomial は β次数24、したがって source scalar c の次数12である。c=1 での多項式恒等式を使ってよいが、これは source 点の評価ではない。

F₄ numerator の B 次数は0、p₁,…,p₄,N₁₂ の numerator の B 次数は1。共通分母へ揃えると、最初の40列の B 次数は高々3、N₁₂² の列だけは正確に4となる。従って QQ 関係があれば最後の係数は0でなければならない。この構造も係数記録に保持するが、rank 判定から N₁₂² の列を除いてはいない。

## 一回の coefficient rank

native FLINT は41個の全整数多項式を生成し、全単項式指数を係数行として集めた。source の点は一つも取らない。係数 union は49,039行。固定素数32003における一回の行消去は

\[
 \boxed{\operatorname{rank}_{\mathbb F_{32003}}=41}
\]

となった。選択した41単項式行の整数41×41行列を全て保存し、その determinant の還元値は **23147≠0**。この minor の整数 determinant は非零なので、QQ でも rank41である。別素数、点の rank、未検証な modular lift を使っていない。

rank40の場合だけ許可された QQ kernel・全係数検証の分岐は実行しなかった。ここで非零 relation は存在しないためである。

## 結論の限界

この41個の universal monomial は独立で、その NEW25 source への制限も単射である。したがって NEW25 方程式が24次なら、この部分空間の外にある。

完全な Inv₂₄ の次元67をこの41個が張るとは仮定しない。従って「全24次式が存在しない」「次数下界26」「NEW25 方程式の決定」は結論しない。追加の部分空間・次数・素数の探索へ自動的に進まない。

## 保存、上限、再現

* `code/q_new_family_degree24_coefficient_relation_20261104_v63.py/.cpp`
* `reports/Q_NEW_FAMILY_DEGREE24_COEFFICIENT_RELATION_20261104_v63/`

`relation.json` は六generatorの順序、全41 exponent、列の定数分母、B 次数上界を保持する。`native-certificate.json.gz` は選択41単項式、全整数41×41係数行列、modular minor、全列項数と検査状況を保存する。full coefficient matrix や全41多項式の重複した保存はせず、hash を固定した六つの exact source coefficient ファイルとコードから再生成できる。

一回の native run、compile と検査込み約16.82秒。最大 polynomial 項数45,994、全272 multiplicationで、60秒・2GB・100k項/多項式・32MiB証明書の上限内。圧縮 native入力約156KB、選択行列証明書約16KB。retry、追加点、追加次数はなし。

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_degree24_coefficient_relation_20261104_v63.py \
  --output work/degree24-coefficient-replay-report \
  --work work/degree24-coefficient-replay-work
```

負の結論は保存した非零係数 minor で再検証できる。actual source の完全な係数・五変数 symbolic section・universal rank6・従来の次数下界の hash を固定し、新しい点での数値的な無関係性とは区別する。
