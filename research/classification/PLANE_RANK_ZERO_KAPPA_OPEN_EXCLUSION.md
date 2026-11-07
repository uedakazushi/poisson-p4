# high-m rank M 二の最後の三源は空：全 normalizer と三つの Jacobi 係数

2026-11-07。主義務 G、producer `plane_collision_review`。
**命題。** 実既約被約 quintic を持つ rank 四 Poisson tensor の rank 零 plane について、
k=1、rank M 二、m≥3 の全 stratum の projective 像は次元≤23。
M(p)非零の四全源の上界は旧 partial packet で閉じ、
その際残した M(p)零 Borel \(\rho=1,5/3,3\) の三つの全源も今回排除する。
rankM三、rankM一以下、k0、他の rank、全 E 又は分類全体は本命題から閉じない。

旧8-file packet は byte を一切変更せず、false の二 report を保持する。
新しい全κopen source の証明は `PLANE_RANK_ZERO_KAPPA_OPEN_NORMALIZATION.md`、
新実行器 `code/classification_plane_rank_zero_kappa_open.py`、
契約 `jobs/plane_rank_zero_kappa_open.json`、
不変証明書 `reports/CLASSIFICATION_PLANE_RANK_ZERO_KAPPA_OPEN_20261107.json`。
gate 後0.484602秒・RSS25,140KiB、全算術 PASS。
本稿は producer 証明であり、独立監査待ちである。

## 1. 全 source と本質的な κopen

旧 first-symbol 全源証明により最後の三源は
\(U=\partial_a\)、\(M=bH+cE\)、\(T=\rho b\partial_a\)、K0、\(\rho=1,5/3,3\)。
全 P1 kernel の係数は順に14,14,15、全 P2 九係数を保つ。
全 tensor は \(\Pi=P_0+P_1+P_2\) である。
旧 flag report の full original350-Jacobi 恒等式は、全 Poisson source 上で

\[
 D|_{w=0}=\kappa z^5,\qquad\kappa=-dr
\]

を与える。d は P1 の \(\{a,w\}\) の z2係数、r は P2 の \(\{b,c\}\) の z2係数である。
κ0 は w が実 D の因子となるので排除済み。
従って全 actual irreducible-D source で d,r は両方非零。
全 parameter の零補集合を失う局所化ではない。

## 2. 実 PGL normalization と有限 source 証明書

前ノートの実変換

\[
 (a',b',c',z',w')=(ra/d^3,b/d,c/d,z,w),\qquad\sigma=d^2/r
\]

は P0 を保存し、d'=r'=1とする。
全150 P0係数の weight零、全 P1/P2 basis の homogeneous weights を直接計算して一致を確認した。
逆は \(a=d^3a'/r,b=db',c=dc'\)、tensor inverse r/d2。
全他係数は自由のままで、全源の forward/inverse は Laurent monomial である。

その全 normalized source に original350-Jacobi を代入すると、
ρ1,5/3 ではまず一つの affine consequence が追加係数ℓ0を与える。
その全350 residual equations から定数1が出る。
ρ3では一段目に定数1が出る。
report は各 current350-row functional、全 parameter map、全350 residual equation を保存し、
元の全350式へ forward map を代入した結果が完全に一致することを検証した。
全 source が空なので、弱い normalized dimension≤10をさらに検査する必要も無い。

polynomial degree は常に二以下、非零 rational pivot だけを用い、
消す unknown の数で制限された affine Gaussian elimination である。
GB、ideal-multiplier degree の探索、jet、Poisson point の抽出は行わない。
旧停止 criterion を繰り返さず、全κopenの実幾何 normalization によって得た矛盾である。

## 3. 全源の短い手計算証明

同じ全係数から三つの Jacobi 係数だけでも全 source の矛盾が分かる。
ここでは P1 kernel の公開順で parameter を p0,p1,… とし、
ρ1,5/3で \(d=p_8,\ell=p_{11},m=p_{13}\)、ρ3で \(d=p_8,m=p_{13}\) と置く。
r はそれぞれ全 parameter 番号20,20,21にある P2係数。
以下の式には全他係数を残してあり、未表示の係数を零にしていない。

全350 row order は quad と degree四 monomial の辞書順で、零始まりの指定は

| row | projective Jacobi quad | monomial |
|---|---|---|
| 51 | (a,b,c,z) | bz3 |
| 121 | (a,b,c,w) | bz3 |
| 122 | (a,b,c,w) | bz2w |

である。\(A=J_{51}\)、\(B=J_{121}\)、\(C=J_{122}\) と書く。

**ρ=1。** 全係数恒等式は

\[
\begin{aligned}
A&=-6r+(2p_0-p_3)\ell+4dm,\\
B&=-3d\ell,\\
C&=-4r+(2p_0-2p_{10})\ell-4dm.
\end{aligned}
\]

d非零と B0からℓ0。
\(A+C=-10r+(4p_0-p_3-2p_{10})\ell\) より r0、矛盾。

**ρ=5/3。** 全係数恒等式は

\[
\begin{aligned}
A&=-8r+(8p_0/3-4p_3/3)\ell+8dm/3,\\
B&=-11d\ell/3,\\
C&=-6r+(8p_0/3-p_3/3-8p_{10}/3)\ell-14dm/3.
\end{aligned}
\]

再びℓ0。
\(7A+4C=-80r+(88p_0/3-32p_3/3-32p_{10}/3)\ell\) より r0、矛盾。

**ρ=3。** 全係数恒等式は

\[
 A=-12r-4dm,\qquad B=0,\qquad C=-10r+6dm.
\]

\(3A+2C=-56r\) より直接 r0、矛盾。
三式とも、全 parameter-polynomial 行列の literal row であり、係数体に generic 条件を課していない。
対応する normalized certificate の定数1 functional は順に
\(-(A+C)/10\)、\(-(7A+4C)/80\)、\(-(3A+2C)/56\) である。
最初の二つでは、その前の B0から得るℓ0の affine map も保持している。

## 4. 全 high-m G への結論と保存範囲

κ0は hyperplane factor、κ非零は §§2–3 の矛盾で、最後の三 source は全て空。
旧 partial proof の common-U rank≤2排除とρ−1の二次 factor を合わせ、
**M(p)零の rankM二 high-m source 全体は空**。
M(p)非零の全四 source は既に≤23である。
従って全 rankM二 high-m stratum が≤23となる。

旧初回の all_Gfalse は「linear kernel bound 単独で全体が閉じなかった」こと、
旧 flag report の false は「定数 rowspace だけではκ零を得なかった」こととして、そのまま真である。
今回それらの bytes や scope を書き換えず、新しい full-source normalization と短い Jacobi 論証を追加した。
全 E 帰属、成分hood・distinctness、新成分数は主張せず、global/main/Git は変更していない。
