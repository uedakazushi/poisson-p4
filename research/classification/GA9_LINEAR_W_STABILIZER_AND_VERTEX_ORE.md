# Ga9 linear-w密軌道：全stabilizerとglobal vertex Ore表示

2026-11-08。PRIMARY E、純手計算。新しいpoint、matrix、profile、jet又はCAS無し。
親成分のpositive membershipはまだ与えない。
canonical tensorと全23軌道閉包について既受理の入力だけから、
許されない退化routeと、全tensorの正しいglobal表示を記録する。

## 1. Connected stabilizerはGa、constant-preconjugate 1PSは使えない

座標(h,x,y,z,w)、h-row零のcanonical四倍tensorは

    Pxy=4h²−y², Pxz=−hx+4w²+2yz,
    Pxw=16hy−4wy, Pyz=hy, Pyw=0, Pzw=−4hw+4y²。

既独立fullambient certificateではEuler24とactual Euler/GL/scalar rank48なので、
projective PGL5軌道は23次元、stabilizer Lie algebraは一次元である。
N=y∂x+h∂zは非零で、各六括弧のLie derivativeを直接取るとL_N P零。
従ってstabilizer Lie algebraは正確にCN。
homogeneous五次元linear matrixのNはN²零、rank二であり、exp(sN)=I+sNは実Ga部分群である。
char零の一次元connected algebraic groupのLie algebraとこのGa inclusionから、
connected stabilizerはこのGaそのものとなる。有限component groupを含めてもGm部分群は存在しない。

任意の固定projective tensor Q、固定g∈PGL5、非自明1PS λ(t)に対し
lim(t→0)λ(t)gQは、存在すればλ全体に固定される。
したがってこのcanonical Ga9点又は同軌道をそのlimitにすることは不可能。
これはT36／T37に限らず、**固定Qと固定preconjugationの全1PS**に対する障害である。

一般のactual arc Q(t)又は非定数preconjugationまで排除したものではない。
親へのmembershipは依然あり得るが、この形のweighted coordinate contractionを繰り返してはならない。
magnetic-zero等の低次元source境界のstabilizerへも外挿しない。

## 2. 同じ点をstrict global vertex Ore表示へ移す

E=h∂h+x∂x+y∂y+z∂z+w∂wとし、同じNを用いて

    Ptilde=P+E∧N

と置く。Euler gaugeなのでprojective tensorもactual Euler-Pfaffianも変わらない。
頂点座標をx、残る四座標を(h,y,z,w)とすると、正確に

    Ptilde=Q+∂x∧W,
    Q=h²∂h∧∂z+2hy∂y∧∂z+(−5hw+4y²)∂z∧∂w,
    W=−hy∂h+(4h²−2y²)∂y+(yz+4w²)∂z+(16hy−5wy)∂w。

QとWはいずれもxに依存しない。Q=∂z∧Vで

    V=−h²∂h−2hy∂y+(4y²−5hw)∂w

はzに依存しないから、Qはstrict Poissonである。
Wの全四成分を使うとL_W Q零となり、したがってPtildeもstrict Poissonである。
vertex rowを単に落としたcoreだけの一致を、元tensorの同値とは呼ばない。
この表示は全raw tensorに対する明示Euler恒等式である。

注意として、先にP+E∧(h∂z)だけを取るとvertex fieldのh成分が零となるが、
そのfieldはQをstrictには保たず、radial defect `E_base∧(−2hy∂z)` が残る。
追加のy∂x Euler項を含めて初めて上のstrict Wとなる。
projective Oreのradial項を無断でstrict zeroと置く欠陥を避けている。

四変数coreのcolumn divergenceは9h∂z、square-zero rank一である。
actual rank四はvertex rowと合わせたものなので、core単独のrank二を
全tensorのrank二と混同しない。
このglobal表示は今後の非定数actual parent arcを組み立てるための正確な入口であり、
現段階で新しいfamily又はknown-component認識を主張しない。
