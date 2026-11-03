# Actual ternary pencil の有限必要locusの独立監査

2026-11-03。`S61_NONSQUARE_TERNARY_PENCIL_FINITE_20261103_v63` の全leading caseを、独自 Fraction 係数演算で検証する。元804項と actual pair の認識は先行独立whole-pair reviewをhash固定で再利用し、著者producerは実行しない。

実ternary二次式の全係数から対称3角行列二つを作り、det(M0+sM3) の四係数を6項determinantとsの係数convolutionで再構成する。保存された四多項式への除算因子がt,P,Jと非零QQ定数だけであることを確認する。元pairの全域rank2は、二式が比例する共通quadratic-factorの場合を既に排除している。残る共通lineがあればpencilの全conicがsingularなので、この四係数の同時零はG2失敗の必要条件である。逆は使わない。

四正規化式の項数28,56,91,131、全次数4,6,8,9を確認し、第一式を `L*t²+M*t+N` とする。全三つのpseudo-divisionを生成し直さず、保存された多項式の積を直接計算し、quotient/remainder恒等式を全係数で照合する。

## L≠0

三式のlinear-t remainderから二つのt消去式F(H,v),G(H,v)を再構成する。取り除いた因子はこのchartの既知unit t,P,J,Lだけである。

保存されたF/GのH係数に対するQQ[v] Bezout和1は、共通vertical factorを排除する。さらにv=0でH次数15,16がそれぞれ保たれることと、QQ[H] specialized Bezout和1を確認する。共通factorがHの正次数を持てば、そのleading coefficientはFとGの両leading coefficientを割る。両leading coefficientがv=0で非零なので、そのfactorのH次数も特殊化で保たれ、specialized和1に反する。従ってF,Gは互いに素で、二変数零点集合は有限である。L≠0では第一式は非零の二次t式なので、各(H,v)上のtも有限となる。

## L=0, M≠0

第一式はM*t+Nである。このlinear polynomialへの保存pseudo-divisionを全係数で確認する。scalar remainderから取り除いた因子はM,JとQQ定数だけであり、このchartで正当である。Lとscalar remainderについて、H係数content和1、v=0でのH次数2,11保持、specialized和1を独立確認する。同じ議論で(H,v)は有限、M≠0によりtは一意となる。

## L=M=0

第一式の零点ではN=0も必要である。L/MとL/Nの小さい固定Sylvester determinantを、独自Fraction permutation和で再計算する。その二つのQQ[v]多項式に対する保存Bezout和が `(v+3)(3v+4)/3` であることを確認する。従ってvは−3または−4/3である。

各vでL,M,Nの保存QQ[H]線形結合が、それぞれH−8、H+14/3の非零定数倍になることを確認する。これでHも強制される。両点で残る三つのactual normalized pencil coefficientをt多項式として評価し、保存QQ[t] cofactorの和が1となることを確認する。このcaseは空である。leading coefficientを全域で逆にしていない。

以上が全場合を尽くすので、aPJ-open 上のsingular-pencil必要locusはQQで零次元以下、その部分集合であるactual ternary共通因子locusも零次元以下である。零次元は空を意味せず、proper finite Z_(p) modelや全点rank78の証明にもならない。撤回済みのG1/G2のみのsyzygy推論は復活させない。訂正版の追加条件G3またはintegral scroll curveを別途扱う必要がある。

再現（新しい出力先）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_ternary_pencil_finite_review_20261103_v63.py \
  --output work/marked-square-v63/ternary-pencil-finite-independent-new
```

全archive/input hashを照合し、60秒・2GB以内で、GB、新しい体、prime、有限体からの昇格、cofactor探索を行わない。
