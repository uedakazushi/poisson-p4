# 被約・連結ファイバーから像の正規性を推論できない反例

2026年11月1日。AB incidence の正規化を調べる際に生じた近道を、明示的反例で除外する。

正規な \(\widetilde B\) から B への proper birational 射の全幾何学的ファイバーが
被約・連結であっても、B の正規性は従わない。
以下の例では始域は滑らか、射は projective、非自明なファイバーは被約な
\(\mathbb P^1\) 一つだけである。

\[
A=\mathbb C[x,y,z]/(xy-z^2),\qquad B=\mathbb C[x,y,xz,yz]\subset A.
\]

A は A1 型の正規曲面錐である。\(z^2=xy\) なので \(A=B[z]\) は B 上有限、
また \(z=(xz)/x\) より分数体は同じ。\(x,y,z\) を次数一とする grading で
\(B_1=\langle x,y\rangle\)、\(A_1=\langle x,y,z\rangle\) だから B は A と異なる。
従って \(\operatorname{Spec}A\to\operatorname{Spec}B\) は非自明な有限正規化である。

B の原点の極大 ideal を \(\mathfrak m_B=(x,y,xz,yz)\) とすると、正規化のファイバーは

\[
A/\mathfrak m_BA=\mathbb C[z]/(z^2)
\]

で、被約でない。

ところが通常の解消
\(X=\operatorname{Tot}\mathcal O_{\mathbb P^1}(-2)\to\operatorname{Spec}A\)
との合成 \(X\to\operatorname{Spec}B\) の全ファイバーは被約・連結である。
二つのアフィン chart を

\[
(x,y,z)=(u,uv^2,uv),\qquad
(x,y,z)=(u'v'^2,u',u'v')
\]

と書く。重なりは \(v'=1/v\)、\(u'=uv^2\)。
原点 ideal の引き戻しは正確に

\[
(u,uv^2,u^2v,u^2v^3)=(u),\qquad
(u'v'^2,u',u'^2v'^3,u'^2v')=(u').
\]

従って合成の原点ファイバーは被約な零切断 \(\mathbb P^1\) である。
原点以外は \(D(x)\cup D(y)\) で覆われ、そこでは z が B の局所化に含まれるため
\(B=A\)。曲面錐の解消もそこでは同型なので、他のファイバーは被約な一点。
解消は projective、正規化は有限なので、合成も projective である。

したがって、正規化のファイバーにある冪零関数が、さらに proper modification を
引き戻すと消えることは実際に起こる。
\(f_*\mathcal O_X\) とファイバーへの base change の単射性を、平坦性や必要な
cohomology の条件を証明せずに使ってはいけない。

AB の proper incidence のファイバーが一点・\(\mathbb P^1\)・LG と被約に分かったことは
有用だが、それだけで AB が正規だとは言えない。
実際の \(f_*\mathcal O\) の生成または局所環を計算する必要がある。

```
.venv/bin/python code/reduced_fibers_normality_counterexample_20261101.py --output /tmp/reduced-fibers-counterexample.json
```

[証明書](../../reports/REDUCED_FIBERS_NORMALITY_COUNTEREXAMPLE_20261101.json)は二 chart の
関係、遷移、原点 ideal の正確な簡約、B の消去による関係式を確認する。
これは AB が非正規という反例ではなく、その正規性に関する一つの推論方法への反例である。
