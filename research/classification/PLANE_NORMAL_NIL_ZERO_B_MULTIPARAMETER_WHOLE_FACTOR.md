# rank零plane・normalnil b零の全mixed／Cartan源：既約quinticでは空

2026-11-08、PRIMARY E。純全source証明と、既存selected-row/kernelの有限検証を区別する。新しいM1、operator、minor、parameter壁又は点を生成しない。後続独立reviewは別receiptで扱う。

## 1. 対象と全normalization

既別作者 `PLANE_RANK_ZERO_M2_COMPLETION_SEPARATE_REVIEW.md` §§1–2 のactual normalizersに従う。normalnil b零、global-vectorが独立な二つの未帰属族は、全higherを保ったまま

```
Π0=T∧(z∂z+w∂w)+V∧z∂w+K∂z∧∂w
```

の次の形へ戻る。

* mixed：H=diag(1,1,−2)、N=b∂a、T=H+τN、V=H+N、K=ηz²、τ≠1。
* Cartan：T=diag(−1−τ,2−τ,2τ−1)/3、V=diag(−1−η,2−η,2η−1)/3、K=κz²、τ≠η。

全係数の零／非零による旧normalizationの補集合は同review及び旧first-jet表へ戻す。このnoteは他のSS型、低parameter二十一源、common-image、scalar又はb非零族を再分類しない。τ=1及びτ=ηはT=Vとなるleading-Pf零のdependent源であり、`PLANE_DEPENDENT_GLOBAL_VECTOR_POINTWISE_COVERAGE.md` の全pointwise E被覆へ戻る。m二の独立族でこれらを無言で削除したのではない。

実Dのleading quadraticは、順に −3b²c(τ−1)z² と −abc(η−τ)z²。従って上記独立条件ではm二であり、位数三以上を同じminorの一般性で扱う必要はない。

## 2. full normal rowの必要条件だけで十分

旧full sourceはweight零35、一30、二9で尽くされる。今回の旧30順序は最初の十八base-pair方向の後に **D,B,E,C** であり、Cはrows27,28,29：

```
D_i z² ∂i∧∂w, B_i zw ∂i∧∂z,
E_i zw ∂i∧∂w, C_i w² ∂i∧∂z.
```

Π0の全z-columnはz倍。weight二P2は全九base-base方向だけで、normal columnを持たない。従って全weight一Jacobiの必要解でC=0なら、実際の全Poisson解でも全z-columnがz倍である。z=0上でraw tensorのimageはこの四次元hyperplane内にあり、Eulerも同じhyperplane内にあるので、signed Euler-Pfaffian Dはzで割れる。Dが零又は可約となり、rank四・被約既約quintic源は空。

ここではP2の特定解、一意性、rank又はhomogeneous kernelを仮定しない。FULL9 P2、全inhomogeneous条件及び全Jacobiを残した必要条件である。

## 3. 保存済みselected kernelの最弱検証

旧producerは各階層点で350×30行列A、r行I／r列J、非零minor、30×(30−r)行列Kを保存している。旧独立reviewはfull operatorと全minor及び壁の階層被覆を検証済みだが、G上界に不要な全kernelのexactnessは一律には主張していない。この違いを保持する。

今回新たに必要なのは次だけである。

1. Jの補集合のK行が単位行列である。
2. 保存されたselected-row行列 A[I,:] K が恒等零である。
3. Kのrows27..29が恒等零である。

旧非零r-minorによりA[I,:]のrankはr、1によりKのrankは30−r。従ってKはselected-row kernelを尽くす。full Aのkernelはこのkernelの部分集合だから、その全元でC=0となる。full A K=0や新rank計算を余分に要求しない。これはfixed slice又は一般点からの推測でなく、当該minor開集合上の全necessary sourceを含むsupersourceの証明である。

## 4. 全壁と二つのsurfaceの閉鎖

immutable四reportの全元階層を、重複を許して次の277件として扱う。

| source | selected-kernel件数 |
| --- | ---: |
| TWOPARAMETERのNIL_B0_MIX | 10 |
| THREEPARAMETERのNIL_B0_CARTAN | 226 |
| HYPERBOLA | 3 |
| CARTAN_SURFACE | 38 |

mixedの旧未閉curve ητ−τ+3=0はη=1−3/τ、τ≠0で尽くされる。τ=0にcurve点は無く、τ=1はdependent。各元linear wallと全有限rootを保持する。

Cartanの旧未閉surfaceは

```
F1=ητ−3η+3κτ+τ,
F2=ητ+2η+3κτ−3κ−2τ−1.
```

actual base swap a↔bとtensor −1は、全higherに作用する可逆変換として (τ,η,κ)↦(1−τ,1−η,−κ) を与え、F1とF2を交換する。F1のτ≠0 graphはκ=−η/3+η/τ−1/3。τ=0ならη=0でT=V、F2の対応するτ=1ならη=1でT=V。従ってvertical κ-linesは既dependent源へ戻る。

CARTAN_SURFACEの九curve graphには、τ=0でcleared行列そのものが零になる保存rootが三件ある。これらは元graphのτ≠0 domain外と明示されており、rank零kernelを実sourceのkernelとして使わない。元surfaceへ戻せば上のτ=η=0 dependent vertical lineだけである。他の274件は全bad-C rowが零。旧独立reviewが確認した全factor／全root／gcd一のgraph分母条件をそのまま使い、新しい壁の抽出又は有限標本への置換をしない。

## 5. 結論と保存範囲

独立T,Vの二whole族では、全parameter、全旧壁、全higherにおいて実Dがz因子を持つ。従って被約既約quinticのE源は空。dependent補集合は既存全pointwise theoremへ戻り、全G上界だけでEを主張していない。

新有限receiptは保存済み277階層のselected-row/kernel必要条件だけを認証する。旧四producer report、旧別作者receipt、埋込producer bytes、旧false及びshared runnerを変更しない。可約Dそのものの全orbit、又は新component数は本命題の対象でない。
