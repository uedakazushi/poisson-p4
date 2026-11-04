# 残るsigmaの実際のfull Poisson射影同値とmoduli体の確定

2026-11-04。候補sigma:(u,v)→(u,−32u−v)は、**全二次Poisson tensorの明示的な有限
射影同値**で実現する。全係数の有理恒等式で確認し、点・GB・高次数不変式は使わない。
前packetの二次の体の挟み込みと合わせると、独立監査後にactual full PGL moduli体の端点が

\[
\boxed{k(X_{\rm NEW25})^{\operatorname{PGL}_5}
       =k(u^2,(v/u+16)^2)}\qquad(k=\mathbb Q,\mathbb C)
\]

と決まる。ここでのclassificationは一般parameterの同定であり、同一parameter内の
有限stabilizer、全特殊parameter、境界の分類ではない。

## 明示的な変換

p=32u+v、epsilon=±1として、元のaffine座標を(x,y,z,w)、新座標をAepsilon(x,y,z,w)+tとする。
Poisson scalarをkappa=−p²/v²とする。

\[
A_\epsilon=
\begin{pmatrix}
0&0&0&8\epsilon v/p\\
-v/p^2&-v^2/p^2&0&8\epsilon v/p^2\\
v(\epsilon v-256)/(512p^3)&v^2(\epsilon v-512)/(512p^3)&-v^3/p^3&4\epsilon v/p^3\\
\epsilon v/(8p)&0&0&0
\end{pmatrix},
\]
\[
t_\epsilon=\begin{pmatrix}
(256-\epsilon v)/(384p)\\0\\
-\bigl(2097153p^3+v(v-512\epsilon)p+2097152v^3\bigr)/(12582912up^3)\\
(p+256\epsilon)/(3072p)
\end{pmatrix}.
\]

det Aepsilon=−v⁷/p⁷であり、表示したopen uv(32u+v)≠0の全体で可逆である。
二つのepsilonのどちらでも、全六bivector成分・全degree0,1,2係数で

\[
\boxed{\kappa A_\epsilon P(u,1,v)(x)A_\epsilon^{\mathsf T}
   =P(u,1,-p)(A_\epsilon x+t_\epsilon)}
\]

がQQ(u,v)上の恒等式として成立する。これは通常のpushforward規約であり、
βのpullback規約と取り違えていない。右辺を新座標の関数へ戻せば
kappa(Ax+t)_*P(u,1,v)=P(u,1,−p)となる。

homogeneous matrix

\[
M_\epsilon=\begin{pmatrix}1&0\\t_\epsilon&A_\epsilon\end{pmatrix}
\]

を使えば、この同じ恒等式はrow/column hを零とした全raw homogeneous degree2 tensorの

\[
\kappa M_\epsilon\,\mathrm{raw}P(u,1,v)(X)M_\epsilon^{\mathsf T}
    =\mathrm{raw}P(u,1,-p)(M_\epsilon X)
\]

という全次数の恒等式へhomogenizeする。従ってEuler quotient後の実際のprojective Poisson
tensorでも等しい。非零scalar kappaはprojectivizationで消える。
σ後のpは−vとなるので、同じopenはσでも保たれる。

## 小さい必要条件からの発見と検査範囲

唯一のtriple H、Q|Hのvertex、既知のunique compatible constant pencilを保つaffine旗を
用い、x'=a x+64c w、w'=c x+a w、y'=e x+D y+g wから始めた。
D=a²−64c²、kappa=1/D、v'=−pとすると必要条件は

\[
j=v'D^2/v,\quad e=D/v-a/v',\quad g=-64c/v'.
\]

二次Pyz比較でyを含む係数は消え、残る三式の組合せは
32ac(u²−256)=0を与える。一般部ではa=0、c=±v/(8p)となる。
残る定数・一次係数は、target modular Gaでy平行移動を零にした後の三未知数についての
五つの線形式で解けた。発見過程のgeneric除外u²−256やv+16uを、最終恒等式の定義域へ
持ち込んでいない。最終変換の検査は上に表示した全openで成立する。

発見receiptはignored work/new-sigma-v63のquadratic.json/full.jsonに保存した。
一次producerはこの探索を繰り返さず、表示した二つの変換をliteral formulaから再構成し、
全六成分を直接検証する。

コンパクト表示への最初の書換えではtzの項−v³/(6up³)を落としたため、Pyz定数係数に
2v³/(9p³)が残って棄却された。この失敗receiptはignored work内に保持する。
上記はその項も含む修正後の式で、元の発見時の厳密解と同じ全係数検査に通過したものだけを採用する。

## 実際のmoduli体の端点とgeneric同定

前packetで、A=u²、b=uv、J=(b/A+16)²に対し

\[
k(A,J)\subset E_{\rm full}\subset k(A,b),\qquad
[k(A,b):k(A,J)]=2
\]

を証明した。known full tensor同値tau:(u,v)→(−u,−v)に加え、今回sigmaも実現したので
全actual PGL invariantはsigmaも保つ。sigmaはAを固定し、b→−32A−bと作用する。
従って

\[
E_{\rm full}\subset k(A,b)^{\langle\sigma\rangle}=k(A,J).
\]

下側包含と合わせて等式となる。これはuniversal invariant環の生成や定義方程式を求めた
結果ではなく、実際のfamilyを通じたfull rational orbit-invariant fieldの同定である。

幾何的generic parameterについて、同じA,Jを持つ可能性は

\[
(u,v),\quad(-u,-v),\quad(u,-32u-v),\quad(-u,32u+v)
\]

の四つだけで、今回のsigmaと旧tauが全てを実現する。したがってこれらが全ての
generic parameter間のPGL同定である。二つのparameter involutionは可換なので、この同定群は
Klein fourとなる。座標変換Mそのものの二乗や、同一tensorの全finite stabilizerを
分類したという意味ではない。特殊uやboundaryに追加同定が生じないとも主張しない。

## 再現と保存

* `code/q_new_family_sigma_finite_equivalence_20261104_v63.py`
* `reports/Q_NEW_FAMILY_SIGMA_FINITE_EQUIVALENCE_20261104_v63/`

証明書はepsilon±双方の全A,t,M,kappa、全六係数の零判定、元full Pと前packetのhashを保持する。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_sigma_finite_equivalence_20261104_v63.py \
  --output work/sigma-finite-replay
```

一回の発見・検証workflowは60秒/2GB内。有限変換の全係数はQQ有理式であり、
新しい点・素数・GB・高次不変式・無制限の非線型equivalence solverは使っていない。
global下側不変量の独立gateと、このfull identityの独立gateが揃ってからmoduli端点を採用する。
