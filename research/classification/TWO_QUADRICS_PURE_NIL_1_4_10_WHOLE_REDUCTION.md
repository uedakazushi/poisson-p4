# H＋相異二既約quadric：純 nil 型一・四・十の全源還元と閉じた部分

2026-11-08、PRIMARY E。既存の全係数証明書を使う純証明。
新しい operator、kernel、rank、minor、profile、jet 又は CAS は実行していない。
新 componenthood／成分数変更は無い。旧全 source、旧 G 証明、主稿、global 台帳、Git は不変。

対象は actual 被約 divisor \(D=hq_1q_2\)、\(q_1,q_2\) 相異なる既約 quadric。
`TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md` の十型を用いる。
以下は型四全空、型一と十の一つずつの閉部分の既知帰属、
各残源を一つの actual orbit へ送る全可逆変換である。
未認識の二軌道について既知親又は新成分を断定しない。

## 1. 型四は H＋二既約quadric を全く持たない

cofactor は \(Z=y\partial_x+z\partial_y+\partial_w\)。
`HYPERPLANE_QUARTIC_NIL31_SHORT_SUPPORT.md` は元全十八自由係数、全九 core 枝、
全 lower／magnetic compatibility を保持する。
`HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_SUPPORT.md` と保存済みの二報告

* `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_SOURCES_20261107.json`
* `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_SUPPORT_20261107.json`

における全十 candidate への入口は、quartic の**既約性ではなく**
最高 degree 四成分が非零という条件だけを用いる。
本対象でも各 \(q_i\) は \(h\) で割れない既約 quadric だから、h一の affine degree は二。
従ってその積の最高 degree 四は非零であり、この全十 candidate 表を適用出来る。

保存 `all_nine_excluded_candidate_certificates` の因子は最高部の因子ではなく、
**全 affine Pfaffian 自体**の普遍因子である。

| zero-based candidate | 保存された全 Pfaffian の因子又は恒等式 |
|---|---|
| 0 | \(z\) |
| 1 | \(z^2\) |
| 2 | \(3cz+1\)、\(c\ne0\) |
| 4 | \(15cz-1\)、\(c\ne0\) |
| 5 | \(z+3t\)、\(t\ne0\) |
| 6 | 最高 quartic 恒等零 |
| 7 | \(z^2\) |
| 8 | \(z\) |
| 9 | \(z^4\) |

最高 degree 四が非零の場合、表の非定数 affine 一次因子は斉次化後にも
actual quartic の一次因子となる。\(3cz+h\)、\(15cz-h\)、\(z+3th\) も全許容 parameter で
非零一次式であり、\(h\) に比例しない。
一方、相異二既約 quadric の積は一次因子を持てない。従って表の九枝は本対象では全て空。
最高 degree が落ちる補集合も candidate 六と同じく単純 \(H\) の条件に反する。

残る candidate 三は \(c\ne0\) の全 \(\mathbb G_m\times\mathbb A^5\) source であり、
同報告が全 parameter に対する primitive な \(x\) 一次式の Gauss 証明を持つ。
その全 Pfaffian は絶対既約 quartic なので、これも二 quadric 積とはならない。

以上で**型四の本対象全源は空**である。
既 H＋既約quartic の MPS 帰属を別因子型へ延長したのではなく、
その前段の全 Pfaffian 恒等式を本対象の正確な因子条件と比較した証明である。
独立入力は `HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_REVIEW.md` と
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_QUARTIC_REVIEW_V2_20261107.json`。

## 2. 型一：全三係数源、T39 部分と一つの残軌道

既 `HYPERPLANE_QUARTIC_NIL_REGULAR_SUPPORT.md` §§2–4、
`TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md` §2 の全 scalar／零補集合を戻す。
\(N=y\partial_x+z\partial_y+w\partial_z\) として、実 \(N^2\)-centralizer 正規化後の全本対象は
\(n_{12}=-1/2,n_8=-1/4,n_{10}=0\)、自由係数
\(A=n_0\ne0,B=n_3,C=n_4\) だけとなる。
元全九十係数へこの代入を戻すと六括弧は正確に

\[
\begin{array}{ll}
 P_{xy}=A+Bz+Cwz-y^2/6,&P_{xz}=Bw+Cw^2-yz/3,\\
 P_{xw}=-wy/3-z^2/12,&P_{yz}=-z^2/4,\\
 P_{yw}=-wz/2,&P_{zw}=-w^2/2.
\end{array}
\]

実 affine Pfaffian は \((2wy-z^2)^2/48-Aw^2/2\)。
\(A\ne0\) なら斉次化は二つの相異なる rank 三 quadric の積となる。
\(A=0\) は同じ quadric の平方として対象外である。

実 \(N^3\)-centralizer \(x'=x-2Cw\)（他座標不変）は、
全六括弧において \(C'=0\) とし、\(A,B\) を変えない。
逆は \(x=x'+2Cw\)。Pfaffian だけの変形ではない。

### B 零：全 T39

既 regular-nil 第二 quartic source は、\(a=n_{10}\) として
\(n_6=3a^2,n_1=-6aC,n_3=0\)、自由 \((a,A,C)\) である。
その六括弧は \(a=0\) にも多項式的で、元 Jacobi は全 \(a\) で恒等零。
\(a\ne0\) 全 source が T39 に入ることは
`HYPERPLANE_QUARTIC_NIL_REGULAR_SUPPORT.md` §§7–8 と
独立 `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_GRAPH_REVIEW_20261107.json` の既受理命題である。
\(a=0\) は正確に上の全 \(B=0\) source なので、閉性から**全 \(B=0\) は T39**。
\(A,C\) の任意の値を保つ actual polynomial source の境界であり、
quartic 既約性を閉補集合へそのまま適用したものではない。

### B 非零：全一軌道への実逆

\(C=0\) にした後、
\[
 (x',y',z',w')=s(t^3x,t^2y,tz,w),\qquad st\ne0
\]
で actual pushforward し、tensor を非零 scalar \(1/t\) 倍する。
六括弧へ直接代入すると同じ canonical 形を保ち、
\[
 A'=s^2t^4A,\qquad B'=st^3B,
\]
かつ \(\operatorname{div}P'=N\) となる。
\(B\ne0\) では
\[
 t^2=A/B^2,\qquad s=1/(t^3B)
\]
を選べば \(A'=B'=1\)。全複素非零値は平方根の二選択で覆う。
逆は各座標を指定された非零尺度で割り、tensor scalar を \(t\) 倍するだけである。
反転した source 係数は既 target unit \(A\) と、この枝の \(B\) のみ。
\(B=0\) は直前で全て処理した。

従って型一の未帰属部分は全 parameter ごとに別の問題ではなく、
**\(A=B=1,C=0\) の一つの actual PGL／scalar 軌道**だけである。
この軌道の親又は componenthood は本稿では未決。
元の \(N^2\)-正規化も明示逆を持つので、元 full source 全体へ戻せる。

## 3. 型十：全三係数源、MPS 部分と一つの残軌道

`HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS.md` §§1–3 と
`TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md` §6 の全 twenty-five coefficient source を用いる。
本因子型は全 affine／完成平方の実逆の後、\(R=x-y^2/2\)、\(\delta a\ne0\)、\(\kappa\) 任意として

\[
\begin{array}{lll}
 P_{xy}=R/2,&P_{xz}=-yz/2,&P_{xw}=-2az+\kappa y,\\
 P_{yz}=-z/2,&P_{yw}=\kappa,&P_{zw}=\delta R
\end{array}
\]

となる。全 Pfaffian は \(\delta R^2/2+az^2\)。
斉次化は \(h x-y^2/2\pm cz h\)、\(c^2=-2a/\delta\)、の相異 rank 三 quadric の積である。
\(\delta=0\) は最高 quartic 零、\(a=0\) は平方として対象外。

### kappa 零：全 MPS

\(\delta,a\) の各非零値を固定し、全実 potential source に
\[
 \Phi_e=\delta R^2/2+z^2(a+ew),\qquad
 V=-R\partial_R/2-z\partial_z/2
\]
を入れる。\(\operatorname{div}V=-1\)、\(V\Phi_e=-\Phi_e\) なので全 Jacobi は恒等零。
元座標の全六括弧は

\[
 (R/2,\ ez^2-yz/2,\ -2az-2ezw,\ -z/2,\ 0,\ \delta R).
\]

全 \(e\) で正則、\(e=0\) は上の全 \(\kappa=0\) source そのものである。
\(e\ne0\) は \(w'=w+a/e\) と既明示尺度により、既認識 translated-minimal C の全軌道へ入る。
同型 C は `HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS.md` §§8–10 の
actual polynomial 曲線から既知 MPS 有限閉和へ入ることが独立確認済みである。
したがって閉性から**全 \(\kappa=0\) source は MPS**。
\(\kappa\ne0\) にこの曲線を使うと \(\kappa e=0\) という Jacobi 条件を失うので、そこへ外挿しない。

### kappa 非零：全一軌道への実逆

\[
 (x',y',z',w')=(t^2x,ty,gz,dw),\qquad tgd\ne0
\]
で actual pushforward し、tensor を \(1/t\) 倍する。
\(R'=t^2R\) から六括弧は再び同じ canonical 形で、
\[
 \delta'=gd\delta/t^3,\qquad a'=td a/g,\qquad\kappa'=d\kappa.
\]
cofactor も \(Z=y\partial_x+\partial_y\) に戻る。
\(\kappa\ne0\) では
\[
 d=1/\kappa,\qquad t^2=\delta a/\kappa^2,\qquad g=td a
\]
により \(\delta'=a'=\kappa'=1\)。
全複素単位を覆い、逆座標・逆 scalar は全て明示的である。
反転した \(\delta,a,\kappa\) の零補集合は上で処理した。
従って残る全 source は **\(\delta=a=\kappa=1\) の一つの actual 軌道**へ入る。
その親又は成分性は本稿では未認識として残す。

## 4. 確認状態と次の最弱義務

元全九十係数、nil 十型、型一の全非線形還元、型四の全十 source／全 Pfaffian、
型十の全 potential map と Jacobi は、各旧 producer と独立 Fraction receipt を入力とする。
本稿の新部分は、それらの literal 全式からの因子条件、全 actual 尺度と逆、
二つの正則全 parameter 曲線による既知閉成分への包含である。
新しい独立算術 receipt 又は新 ambient computation を実行したとは主張しない。

型一の \(B\ne0\) と型十の \(\kappa\ne0\) は、二つの finite representative へ戻っただけで未被覆である。
もし既存 faithful S61 chart／既知閉 source への陽性認識を使うなら、
各代表が実際に入る symplectic chart、元全 tensor への正逆、全 graph 方程式を必要とする。
一つの chart での陰性を閉成分非所属と呼んではならない。
型七 full twenty source、他の factor、全 E、総分類は別に保持する。
