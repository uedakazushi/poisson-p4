# reduced R・m四の残六全源：実 shear 全逆と因子の独立確認

2026-11-08、PRIMARY E。**独立 PASS：指定六 fixed 全源の被約既約 quintic locus は空。**
旧 G 次元上界を E 帰属へ読み替えず、全 source 上の実正規化を確認する。
新成分数は零。他の m四分岐をこの六つだけから尽くしたとは主張しない。

## 1. 全 source の適用範囲

既独立 special-j／repeated-leading review の正常形・actual逆・全必要核を利用する。
対象は `J0_secant,J0_compression,I0_mixed,triangle_diagonal,triangle_Jordan,tangent_conic_net`。
前五は triangle、最後は tangent conic+line の fixed first-symbol。
座標は `(a,b,c,u,v)`、P0は元の全150係数をそのまま使う。

raw weight一は36方向、Euler六方向を除く完全源は三十方向。
独立 recipe で全30 basisと全9P2 basisを作り、旧保存全150 coefficient vectorsと照合した。
P2は base pair ab/ac/bc ごとの u²,uv,v²、全九方向であり normal column は無い。
旧 K1 は全 Jacobi weight一＋actualD3零に対する complete 六核である。
一般元の軌道次元だけから shear transitivity を推測していない。

## 2. 全 six-shear 正逆を全十五 parameter で確認

actual map は

```
a'=a+s0u+s1v, b'=b+s2u+s3v, c'=c+s4u+s5v,
u'=u,v'=v。
```

det一、inverseは全sの符号反転。
weight一 Euler商の独立な全域 inverse も直接用いた。
raw P1 の `P_i,u[u²]` と `P_i,v[v²]` を読み、
`Y_i=-P_i,u[u²]*u-P_i,v[v²]*v` とすると、P1-E∧Yが指定30 basisへ一意に戻る。
残る30係数の再構成が元全 raw tensor と一致することを確認した。
新しい quotient rank 又は kernel を解いていない。

六つの unit shear で元 P0 の全150係数を直接 pushし、保存 S(30×6) と一致。
旧 K(30×6) と全 constant C について `S=KC`、Cと保存 inverse の両積は単位行列。
全六 source でこの等号が成り立つ。

さらに unit six 列だけで止めず、**全α0,…,α5と rawP2九係数を同時に自由変数**として

```
P=P0+P1(Kα)+P2(r), s=-C^-1 α
```

の全 actual pushforward `G*P(G^-1x)*G^t` を別 Fraction で計算した。
全 parameter 上で weight零は同じ P0、weight一は Euler quotientで零。
残 weight二は正確に全9P2へ入り、その係数は

```
r_new=r+q(α)
```

という global triangular map である。全九 q(α) を report の疎多項式で保存した。
qは元 r に依存せず、inverse `r=r_new-q(α)` は全域 polynomial。
これにより affine inhomogeneous particular を含む全 rawP2 を保持することを、
実際の全 source 恒等式で確認した。source parameter の反転は一つも無い。

normal u,vと plane上のbase座標は shearで固定されるので、
normal order四の leading shape 条件も変わらない。
base substitution が生む差は更に高い normal orderであり、G4自体は同じである。
従って正規化された全対象の r_new は旧全 weight二/leading-shape matrix L の kernelへ入る。

## 3. P1除去後の全9P2・全350・actual D

P1を**実際に**消したため二次 forcing が零となり、正しい条件は L*r_new=0。
homogeneous kernel零を P2そのもの零と仮定する議論ではない。
保存 normalized9×n kernelを元 L の全363／364行へ戻して零を確認し、
旧 rank certificate の自由座標上で単位行列になることを確認した。
旧の正確な rank と合わせ complete kernel であり、新しい kernel 計算は行っていない。

その全 kernel係数 t_i を自由に残した P0+P2について、
全十括弧／150係数、元全350 Jacobi、定義からの full signed Euler-Pfaffianを別に再構成した。
表示 source と全係数で一致する。
actual Dは normal degree四のみで、normal degree五の未処理項は無い。

| fixed 全源 | 完全P2核次元 | 実Dの十分な因子又は零 |
|---|---:|---|
| J0_secant | 0 | 恒等零 |
| J0_compression | 2 | uv |
| I0_mixed | 1 | uv³ |
| triangle_diagonal | 1 | u²v² |
| triangle_Jordan | 2 | u²v |
| tangent_conic_net | 1 | u² |

元 D と保存された全 factor quotient の積も正確に一致。
全t_i零／非零を含み、非零Dはproper positive-degree factorを持つ。
u,vは全actual shearとそのinverseで固定されるので、元 sourceへ戻っても同じ因子が残る。
可約Dのrank四点までrank二と呼んだものではない。

## 4. 検証と残る scope

ONE finite gate 後、Fractionのみ **0.174秒／RSS28560KiB／206恒等式**。
producer/CAS import、新M1、rank、kernel、point、profile、GB無し。
全 simultaneous15parameter pushと P2 shiftsを保存しており、単位shearの線形性だけの外挿ではない。

- job：`jobs/plane_m4_six_whole_separate_review.json`
- receipt：`reports/CLASSIFICATION_PLANE_M4_SIX_WHOLE_SEPARATE_REVIEW_20261108.json`

既受理の六正常形完全性と併せ、reduced R・m四のこの残六 E 義務は全て空となる。
他の rank一、partition22/4、concurrent、integral R、nonreduced R の全被覆は別の受理済み入力又は別義務。
全 rank二plane m四の aggregate と台帳登録は統合担当がそれらを合わせて行う。
producer／旧source／currentjob／主稿／global／Gitを変更していない。
