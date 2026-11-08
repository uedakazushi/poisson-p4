# reduced R・m四：残六fixed接続の全実源は既約quinticで空

2026-11-08、`/root/quartic_review`。PRIMARY E。
既受理の全source幾何・全operator・完全kernelを保持したまま、残六源について
**全P1を実shearで消せる**ことを証明し、全P2と実Pfaffianから空性を得る。
新しいfirst symbol、M1、点、profile、normal jet又はGBは用いていない。
本稿の新算術・全源逆の別作者reviewは作成時点では未了。

## 1. 対象と全補集合

対象は既root-cover/pencil正常形の次の六つだけである。

| R/leading | report type |
|---|---|
| triangle・SF J0 secant | J0_secant |
| triangle・SF J0 compression | J0_compression |
| triangle・SF I0 mixed | I0_mixed |
| triangle・211 semisimple pencil | triangle_diagonal |
| triangle・211 Jordan pencil | triangle_Jordan |
| tangent conic+line・211 net | tangent_conic_net |

squarefree三源は既 `PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_SEPARATE_REVIEW.md`、
211三源は既 `PLANE_M4_REPEATED_SEPARATE_REVIEW.md` の完全kernelを用いる。
concurrent既空源、他のm四rank一／22／4／非被約Rの既知八成分被覆を再開しない。
元の仮定はactual D既約被約quintic、restriction rank二plane、normal位数四。
全leading amplitudeの零を大きな必要源へ残すので、それらを勝手に単位化しない。

## 2. 六shearと完全K1の全域同型

全座標は(a,b,c,u,v)。各fixed P0に対し全higher源をP0+P1+P2と書く。
P1はraw36/Euler6の完全30基底、P2はbase pair三つ×(u²,uv,v²)の**全9**。
元Jacobi weight一とactualD3零の完全kernel K1は六源ともdimension六である。

実座標変更

    a'=a+s0u+s1v, b'=b+s2u+s3v, c'=c+s4u+s5v,
    u'=u, v'=v

は全s上でdet1、逆は全sの符号反転。
full quadratic tensorのpush後、weight一Euler6を一意に除く。
raw36上の全30基底＋Euler6の定数逆は既疎Fraction基底から作る。
六個のunit shearの全150pushと、六s同時のpolynomial pushを計算し、
weight一作用が保存した30×6行列Sの厳密な線形作用であることを確認した。

旧完全kernelの列をK（30×6）とし、P1=Kαと書く。
全六源で **S=KC、det C≠0** で、全定数逆α→sを保存した。
最初の五源では、順にη=8/3,4/3,2,4/3,−2/3として

    (s0,s1,s2,s3,s4,s5)
      =η(α3,−α0,α4,−α1,−α5,α2).

tangent-netでは

    s0=4α0/9, s1=4(α3−4α1)/9, s2=4α1/3,
    s3=4(2α2+α4)/3, s4=−4α2/3, s5=−4α5/3.

αの順は各旧full30_weight1_kernelの保存順である。
これはdim K1=6からtransitivityを推測したものではなく、全係数の実線形同型である。
どのsource parameterも反転せず、全零補集合にそのまま適用出来る。
逆shearを全tensorへ実施すると、P1は零、P0は同じfixed P0である。
変換から生じるweight二項は全9P2に入る。
したがってinhomogeneous P2 particularを落としていない。

## 3. 全9P2と全actual D

逆shear後はP1=0なので、元weight二の全Jacobi＋全leading shape条件は

    L(P2)=0

となる。これは元inhomogeneous方程式を無根拠にhomogeneous kernelへ狭めたものではない。
**実際にP1を消したため**quadratic forcingが零なのである。
全9未知量への既完全Lのkernelを使い、全normal-root振幅を自由に残す。
P2はbase-baseのnormal二次式だけなのでJ(P2)=0、P2∧P2=0。
各normalized全源の元350projective Jacobiも直接零と確認した。
actual PfaffianはP0とP2のcrossだけ、normal degree四でありg5項は無い。

全kernel係数をt_iとした結果は次の通り。

| type | dim ker L | normalized actual D |
|---|---:|---|
| J0_secant | 0 | 0 |
| J0_compression | 2 | `−3uv(−a t1u²+b t0v²)/4` |
| I0_mixed | 1 | `−b t0uv³/2` |
| triangle_diagonal | 1 | `−3b t0u²v²/4` |
| triangle_Jordan | 2 | `3u²v(−a t1v+b t0u)/2` |
| tangent_conic_net | 1 | `9t0u²(a u²+2buv+c v²)/4` |

全t_iは任意で、非零／零を捨てていない。
第一源の全actual Dは零、他は正次数normal factorを持つか零。
normal u,vは六shearで固定されるため、逆座標変更を戻した元Dにも同じproper factorがある。
したがって全六源の既約quintic locusは空。
可約Dを持つrank四点までrank二と呼ばない。全rank四空と証明したものでもない。

## 4. 算術receiptと既結果への接続

有限契約 `jobs/plane_m4_six_whole_shear_factors.json` は実行前gate PASS。
exact実行programは `inputs/plane_m4_six_whole_shear_factors_20261108.json` に保存し、
新shared appendやwrapper codeは作っていない。
旧focused sparse5/constant shear engineを不変のまま再使用した。
定数逆とpolynomial係数にはSymPyのexact QQ演算を用い、数値rankや有限体は用いない。
全六一回jobは4.046998秒、RSS64,944KiB。

report `CLASSIFICATION_PLANE_M4_SIX_WHOLE_SHEAR_FACTORS_20261108.json` は
全30×6S、全CとC逆／det、kernel座標row、全9P2 kernel、normalized全十括弧、
元350残差、全actual Dと因子商、実行program及び全input hashを保存する。
旧M1／K1／source／report／reviewerのbytesは改変していない。

既幾何の六正常形尽くしと合わせれば、reduced Rのm四残六E義務は全て空になる。
既rank一／22／4等の八成分被覆を合わせたaggregateは統合担当が行う。
nonred Rの未割当義務やm二／三／全研究の完了へ本稿だけを拡張しない。
