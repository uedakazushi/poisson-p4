# CA V1の数値Lie HessianとQQ一変数source jets

既存のcheap checkpointを再利用し、symbolic多変数Hessianを作らずに、全150×29の \(H_{\kappa,j}\) を得た。実行時間は0.214秒で、全29列の二次polarizationとdiagonal作用の7次までのclosed weight式に一致した。57 block、八つのactual odd制約、その後のcusp orderはこのhelperの結論ではない。

## sourceとgroup chart

source座標順は、checkpointに保存された20個のGL generator、次いで九個の固定fiber方向である。

\[
F(s,t)=\rho\bigl(\exp(s_0E_0)\cdots\exp(s_{19}E_{19})\bigr)
       \bigl(P+\sum_{a=0}^8t_aF_a\bigr).
\]

各 \(E_i\) は保存lex順のelementary matrixである。\(\rho\) はquadratic bivectorの実際のpushforward作用で、150係数の整数疎行列 \(L_i\) で微分できる。出力tensorの二つのindexへの作用と、二次monomialの変数置換微分を係数ごとに構成した。

以前の失敗したrunnerは \(I+s_iE_i\) を使っていた。今回のchartは **formal exponential** である。offdiagonalでは一致するがdiagonalの二次以上は異なる。一次Jとkernelは同じでも、古いdiagonal Hessianを流用しない。

## ordered Hessian

このgroup順序に対し、group二変数の二階微分は

\[
H_{ij}=L_{\min(i,j)}L_{\max(i,j)}P
\]

であり、diagonalも \(L_i^2P\) である。group/fiber mixed微分は \(L_iF_a\)、fiber同士は零である。これを保存kernel \(\kappa\) とcontractして全29列を得た。150係数のJも係数作用から作り直し、checkpointと厳密に一致し、\(J\kappa=0\) であることを確認した。

保存先は `reports/AMBIENT_CA_V1_COEFFICIENT_HESSIAN_20261103_v63/hessian.json.gz` で、主要keyは

- `H_kappa_original150_by29`
- `J_original150_by29`
- `kernel29`
- `complement28_columns`
- `kernel_kernel_quadratic150`

である。全て元のsource29順序であり、使用側で \([\kappa,S_{28}]\) へ変換する。

## 一変数jet evaluator

同じmoduleの `evaluator(point150,basis9_columns,group_GL_indices,degree)` は、29本のQQ係数source seriesを150本の出力seriesへ送る。source seriesの定数項は零である。

group積の作用は右端のfactorから適用し、各factorについて \(\sum s_i^kL_i^k/k!\) を必要次数まで厳密にconvolutionする。\(s_i(0)=0\) なので有限項で足り、positiveな最低次数に応じて不要な項を省く。diagonalのformal exponentialを切り捨てたり、数値浮動小数点に置き換えたりしない。

二次では全29個の \(\kappa+e_j\) をpolarizeし、各Hessian列を検査した。\(\kappa\) 一方向だけの確認ではない。また20 generator中の四diagonalについて、各150 monomialのweight

\[
\mathbf1_{u=a}+\mathbf1_{v=a}-\#\{a\text{ in monomial}\}
\]

によるclosed exponential式と、非線形source series \(s=\epsilon+\epsilon^2/3\) の次数7までを全係数で比較した。

これらは係数Lie表示内での厳密検査である。八actual odd制約や57 blockを解く側は、別の5×5有限座標source evaluatorとの比較も行い、その後の形式陰関数補正を全150出力で検証する必要がある。helperだけからcusp型を宣言しない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_CA_V1_coefficient_hessian_20261103_v63.py \
  --checkpoint reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz \
  --output work/CA-V1-coefficient-Hessian-replay
```
