# 階数三の混合余剰クラス：新しい22次元成分の独立証明

2026-11-07。主課題 E_quadric3。全射影 Poisson 空間で、以下の実際の
PGL5 軌道閉包が22次元既約成分であり、従前の58成分と相異なることを確認した。
したがって成分数の下界は59となる。全空間の網羅性は依然として未証明である。

## 1. 点と独立な全350式

座標を `(h,u,v,z,w)` とし、h を Casimir にした二次斉次二ベクトル P を
上三角の六係数の順 `(uv,uz,uw,vz,vw,zw)` で

\[
 (-hv,\ 3hz-2uw,\ 2hw+2u^2,\
 2hw+2u^2+2wv,\ -2uv,\ 4w^2+6uz)
\]

と定める。raw Jacobi の全十成分は零で、

\[
 \operatorname{Pf}(P)=4(hw+u^2)^2,\qquad
 D_P=4h(hw+u^2)^2.                              \tag{1}
\]

Q=hw+u² の階数は三、従って Q は既約であり、P は一般階数4である。

独立実装は、全150個の二次斉次係数を動かし、全350式

\[
 F(P)=E\wedge\operatorname{Jac}(P)=0
\]

の微分 D を整数で再構成した。`Jac` は通常の三項循環和とする。
固定 h-Casimir 部分や固定退化因子の接空間ではない。
producer の係数演算・単項式列挙・Jacobi・群作用・階数 helper は一切 import していない。
単項式は指数ベクトルを独立に列挙して逆辞書式順に並べた。

V を全25列の Euler-vertical 方向、M を V・全25 GL 微分・scalar 方向の結合行列とする。
全て実際の群作用・線形商の方向である。次の R を追加した。

\[
\begin{array}{lll}
 R_{hu}=-2hu+2uv,&R_{hv}=-4hv+6v^2,&R_{hz}=4hz-2vz,\\
 R_{hw}=2u^2,&R_{uv}=-uv,&R_{uz}=3uz,\qquad R_{vz}=11vz,
\end{array}
\]

他の上三角係数は零である。全行列を QQ で消去し、次を独立に確認した。

| 行列 | QQ 階数 | 非零小行列式 |
|---|---:|---|
| 全350×150微分 D | 102 | `3756079923488420961319000576855418466046094478409728000000000000` |
| V | 24 | `1` |
| M | 47 | `18874368` |
| `[M R]` | 48 | `−75497472` |

小行列の全行列番号は report に保存した。さらに **`D[M R]=0` を整数で検査**した。
従って R は群作用・鉛直・スカラー方向の外にあり、全接空間の残る一方向を張る。
射影軌道次元は `47−24−1=22`、全射影接空間次元は `150−102−24−1=23` である。
接空間の余剰一方向だけでは新成分を結論せず、次の障害を使う。

## 2. 十項の余核と真の一変数形式切片

全350方程式上の整数線形汎関数 λ を次の十係数で与える。
添字は `(h,u,v,z,w)=(0,1,2,3,4)` に対応する。

| 四ベクトル成分 | 単項式 | λの係数 |
|---|---|---:|
| 0123 | h²uz | −18 |
| 0123 | h²w² | 18 |
| 0123 | hu²w | −9 |
| 0123 | huvz | 9 |
| 0123 | uv²z | −5 |
| 0124 | h³z | 36 |
| 0124 | h²uw | 18 |
| 0124 | h²vz | −36 |
| 0124 | huvw | −9 |
| 0124 | uv²w | 4 |

独立に

\[
 \lambda D=0,\qquad \lambda F(R)=-140\ne0       \tag{2}
\]

を確認した。従って任意の二次補正 S を許しても

\[
 F(P+tR+t^2S)=t^2\bigl(F(R)+D(S)\bigr)+O(t^3)
\]

の t² 項は消せない。方程式係数の並びは上表で物理的な四ベクトルと単項式へ戻しており、
producer の行番号を意味づけずに流用したものではない。

47個の実際の群・鉛直・scalar 方向に横断する103次元の線形切片を、R を接方向に含めて選ぶ。
全 D の102次可逆小行列式による形式陰関数定理を適用すると、真の残余環は

\[
 A=\mathbb C[[t]]/J
\]

となる。残る248方程式を捨てていない。
その全方程式のうち少なくとも一つは、(2)によって非零二次初期項を持つ。
従って J≠0 であり、`dim A=0`。

群作用と横断切片による形式座標、または47個の切片方程式に対する次元不等式を戻すと、
全150係数空間での局所次元は高々47となる。
一方、実際の群・鉛直・scalar 軌道は47次元である。
従って局所次元は47、鉛直商と射影化の後では22となる。
22次元 PGL 軌道の閉包を含む全 Poisson 成分の次元も22以下なので、
この軌道閉包自身が全空間の22次元既約成分である。

さらに P を通る任意の既約成分は、連結 PGL5 が各成分を保つため、その22次元軌道を含む。
局所上界22から、全てこの軌道閉包に一致する。
従って P を通る被約成分は一つである。
同じ点を含む既約な実際の族の像閉包は、この成分に含まれる。
これは全 mixed 族の支持を接続する際に使えるが、任意の unrelated 族への主張ではない。

ここでは被約支持と成分性に必要な零次元性だけを使用する。
正確な非被約厚み・一般横断長・正規性を追加の結論とはしない。

## 3. 従前58成分との相違

MPS40 の射影次元は、既存独立 crosswalk で23が29件、24が8件、25が3件である。
従って新しい22次元成分と相異なる。
追加11成分と NEW25 のうち、22次元である既知一般 rank4 成分は
`DOUBLE_SMOOTH_QUADRIC` と `TRIPLE_HYPERPLANE_QUADRIC` だけで、
他は全て23次元以上である。二つの22次元成分とは以下の閉条件で分離する。

* 射影線形形式 L と、階数三以下の二次形式 Q から `LQ²` を作る掛け算写像は、
  射影的な源からの射なので、その像は閉である。零を含む affine cone を取れば、
  退化五次式が零になる rank2 境界も保持できる。
  (1)の全 PGL 軌道とその閉包はこの条件を満たす。
  既知 `DOUBLE_SMOOTH_QUADRIC` の一般元は `Lq²`、rank(q)=5であり、
  既約因子分解の一意性から階数三以下の二次因子へ書き換えられない。
  従ってこの既知成分とは一致しない。
* 「三乗の一次因子を持つ五次式」は同様の射影的掛け算から閉条件となる。
  既知 `TRIPLE_HYPERPLANE_QUADRIC` の全閉包はこの条件を満たす。
  (1)は、相異なる既約因子 h と Q の重複度が1と2なので、三乗の一次因子を持たない。
  従って P はこの既知成分に属さない。

既知六 rank2 成分の全閉包は `π∧π=0` という閉条件に入るが、(1)は非零である。
従って新成分はそれらとも相異なる。以上により rank4 の下界は53、全体の下界は59になる。
次元表の元定理・旧証明書は、既存の早期・後期成分監査と NEW25 監査への依存である。
今回新しく再構成したのは、この新点の全350式、全群作用と障害である。

## 4. 証拠と fresh clone での再現

* 独立実装：`code/classification_rank3_mixed_component_review.py`
* 独立 report：`reports/CLASSIFICATION_RANK3_MIXED_COMPONENT_REVIEW_20261107.json`
* producer の永続保存：`reports/CLASSIFICATION_RANK3_MIXED_COMPONENT_20261107.json`
* 契約：`research/classification/jobs/rank3_mixed_component_review.json`

独立実行は約1.75秒、制限60秒・1.5GiBで PASS。
producer の永続 JSON は、独立実行時の `work/rank3-mixed-support-obstruction.json`
と **同一 bytes** で保存してあり、SHA256 は
`0e887e60e8f78fccd7f575d9e7c9def491de1be1d1e10ebc69bbd8c0e635cb69` である。
旧 review の数値・code hash・入力 hash を変更せず、元入力パスを次で復元できる。
producer JSON 内には全 source code も保存されている。

```sh
mkdir -p work
cp reports/CLASSIFICATION_RANK3_MIXED_COMPONENT_20261107.json work/rank3-mixed-support-obstruction.json
python3 scripts/classification.py gate research/classification/jobs/rank3_mixed_component_review.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python code/classification_rank3_mixed_component_review.py --output work/rank3-mixed-independent-replay.json
```

出力は既存ファイルを上書きしない。
独立実装は producer を実行せず、literal な P と R から全行列を作り直し、
保存した十係数の余核も物理的係数の順序で検証する。
全局所環・全未知成分・他の退化因子の支持を尽くしたという主張ではない。
