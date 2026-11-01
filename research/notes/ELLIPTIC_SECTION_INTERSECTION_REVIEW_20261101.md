# 二つの楕円グラフの交叉：独立点検

2026-11-01。T39、T40 の被約閉包を既証明のグラフ Γ₅₁、Γ₅₂ とする。
**そのスキーム論的交叉は、5次元・次数24の rank-one flag cone を被約台とし、
頂点に長さ10の二乗零な nilradical を持つ。**
これは二つの被約成分の交叉の記述であり、周囲の全S61局所環を分類したものではない。

## 1. 正確なイデアルの意味

R=QQ[W16]、m=Rの正次数イデアルとし、

    I=(Q₅₁−Q₅₂),
    J=(primitive β の5×4係数行列の2×2小行列式)

と置く。二つのグラフの定義イデアルの和を取り q を消せば、交叉の座標環は正確にR/I。
この消去には飽和、一般点への制限、根基化を必要としない。

独立コードは新しい交叉計算コードを import せず、保存された二つのsectionと既存の
primitive β 構成から I、J を作り直した。元の検証とは逆の最小monomial pivotで有理
線形代数を行い、次を確認した：

    dim I₂=45,   dim J₂=55,   I⊂J,
    dim I₃=560,  mJ⊂I.

両イデアルが二次生成なので、ここでの次数2・3の確認は記したイデアル包含全体を証明する。

## 2. Hilbert級数の別方式による全次数確認

外部CASから得た55個のGröbner候補を、検証されていない答えとして入力した。
その次数2の張る空間が J₂ と一致することを確認し、互いに素な先頭項を持つ対には
product criterionを適用、残る376個のS多項式を全て再還元して零を確認した。

元の検証はcolon idealによるHilbert再帰を使うが、本点検は別の数え方を使った。
先頭二次単項式は、53個の異なる変数の積と二つの平方である。
16変数の全65536個のsupportを列挙すると、許容されるsupportは348個。
平方で制約される二変数の指数は0か1で、他のsupport変数の指数は任意の正整数なので、
supportの寄与は

    t^(free+bounded)/(1−t)^free

となる。これを全supportで足すと、有限個のHilbert値の照合ではなく有理関数の恒等式として

    H_R/J(t)=(1+11t+11t²+t³)/(1−t)⁵

を得る。

## 3. rank-one locus の幾何と J の素性

零でないrank-one tensorは β=ℓ⊗ω と書ける。
閉性 ℓ∧ω=0 は ω=ℓ∧m を意味し、primitive条件は二平面 span(ℓ,m) のisotropic性を
意味する。よって射影化された被約台は flag variety F(1,2;4)=Sp₄/B。
埋込みのhighest weightはω₁+ω₂である。

このhighest-weight coneの次数n座標空間はCartan power V(n,n)*。
Weyl次元公式によりその次元は(n+1)⁴で、前節と同じHilbert級数になる。
座標環のCartan-power記述は、一次資料
[Johnson, Multivariate Hilbert series of lattice cones and homogeneous varieties](https://bpb-us-w2.wpmucdn.com/sites.uwm.edu/dist/d/129/files/2016/04/wayne-a-johnson-1wx801b.pdf),
§5.3.2 Theorem 5.3.7の一つのweightの場合にも記されている。

Jはflag coneの素イデアルPに含まれる。従ってR/J→R/Pは次数を保つ全射である。
両辺の**全Hilbert級数**が等しいため、その核は零でありJ=P。
単に次元5、次数24が一致するだけで素性・根基性を結論したのではない。

## 4. 頂点における二乗零部分

I⊂J、mJ⊂Iと二次生成性により、J/Iは次数2だけにあり、

    J/I ≅ QQ(−2)^10

というR加群になる。J⊂m²だから J²⊂I。従ってJ/IはR/Iの二乗零イデアルである。
Jが素なので、これが正確にR/Iのnilradicalとなり、被約化はR/J。
nilradicalはmに消されるため頂点にのみ支持され、長さは10。
頂点以外では交叉は被約で、flag coneの滑らかな開部分に等しい。

ここでの加群の同型は、R/IがR/Jと二乗零空間の直積として分裂するという意味ではない。
「追加の10点」という集合論的解釈もしない。
関連素イデアルは最小素イデアルJと埋込み素イデアルmだが、全S61の関連素イデアルを
列挙したものではない。

```sh
.venv/bin/python code/elliptic_section_intersection_review_20261101.py \
  --output /tmp/elliptic-section-intersection-review.json
```

証明書：`reports/ELLIPTIC_SECTION_INTERSECTION_REVIEW_20261101.json`。
逆pivotによる包含、376対の再還元、348supportの個数表、別方式の全Hilbert級数を保存した。
