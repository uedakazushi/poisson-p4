# dependent normalnil・trace零全源は二つのresonance以外G空

2026-11-07、主義務G。producer手証明、独立review前。
scalar/traceless零の別noteではなく、ここではglobal V≠0を保つ。

全重複度のmodular入口から、rank零planeのfull35/30/9 sourceを使う。
dependent pairのnormalnilでtrace field T=0の場合、全weight零は

`Π0=zV∧∂w+αz²∂z∧∂w`、`V(x)=A*x`、trA=0、A≠0。

全Π1三十、Π2九を残す。
Π1のnormal方向の一項を`w² C∧∂z`、C∈C³ constant-base vectorとする。
他の係数を零にしたsliceを選ぶのではなく、全Jacobiの一つのnormal monomial coefficientだけを見る。

## 1. 全JacobiからCのnecessary identity

projective Jacobiのbase二・normal二quad成分で、normal monomial z*w²の係数は正確に

`x∧(AC+2αC)−3(Ax)∧C`。                              (1)

この係数へ他のΠ1項又はΠ2は寄与しない。
実際zL,wM、Bzw、Ezwの寄与はz³又はz²w、D z²の同係数は零。
full higher Jacobiはweight二以上なので、このweight一係数を変更しない。

(1)のV部分は直接、strict bracketの

`[Π0_V,Π1_C]^{ijz}=2zw(V∧C)ij`、

`[Π0_V,Π1_C]^{ijw}=−w²(V∧C)ij`、

`[Π0_V,Π1_C]^{izw}=zw²(AC)i`

をEuler wedgeへ戻して得られる。
K部分のstrict izw係数は2αzw²Ciなので、追加は2αx∧C。
これはcanonical sourceの全係数比較であり、global divergence liftを選ぶ議論ではない。

## 2. C零ならactual Dはlinear factorを持つ

全canonical sourceでΠ2はbase-baseだけ、vertical normal係数はαz²である。
従ってz=0上でnormal contractionπ#dzに残る項は符号を除きw²Cだけ。
C=0ならhyperplane H=(z=0)はprojective Poisson hypersurfaceとなる。
Hのtangent dimensionは三なのでπのrankはH上で二以下、actual Pfaffianはzで割れる。
Gのactual Dは被約既約quinticなので、この場合は不可能。
ゆえにGではC≠0が必要である。

## 3. C非零の場合の全matrix代数

(1)を全xで零としx=Cを代入すると4C∧AC=0、従ってAC=μC。
再び全xに戻すと`((μ+2α)I−3A)x`のimageはC-lineに入る。
よってあるlinear form ℓについて

`A=((μ+2α)/3)I+C⊗ℓ`。

trA=0からℓ(C)=−(μ+2α)。
AC=μCも使うと`μ=−2(μ+2α)/3`、したがって

`μ=−4α/5`、`A=(2α/5)I+C⊗ℓ`、`ℓ(C)=−6α/5`。          (2)

α非零ならAはC-line上の固有値−4α/5と、quotient上の固有値2α/5を持つ。
両者は異なるのでAはsemisimple double-eigenvalue型であり、mixed Jordan型ではない。
α零ならA=C⊗ℓ、ℓ(C)=0なので、A≠0からrank一nilpotent型だけである。

## 4. 全Jordan表への帰結と正確な残余

`PLANE_RANK_ZERO_DEPENDENT_GLOBAL_VECTOR_FIRST_JETS.md`のnormalnil trace零全源に適用する。
regular semisimple、mixed、regular nilpotent Vは全αでG空。
rank一nil V=Nではα非零でG空。
semisimple double V=H=diag(1,1,−2)では、quotient固有値1と(2)を比較して
`α=5/2`だけが残り、Cは固有値−2のlineに入る。

したがって今回の全T零normalnil sourceで未処理として残すのは次の二whole sourceだけ。

- `V=H`, `K=(5/2)z²`、Cはdistinguished eigenline、C非零。
- `V=N`, `K=0`、Cはimage N、C非零。

この必要条件だけで両残余の全Poisson support又はG非優越を結論しない。
full thirty/nine operator又は全original higher Jacobiを使う有限課題として保持する。
前者のaffine first-symbol base上界は六（trace零H源のbase七からα固定一）、後者は五。

T非零のdependent源、normalssのtrace零、scalar V零、又は全高mの完了には拡張しない。
既知低rank familyへ退化するだけで元sourceのmembershipを推論することもしていない。
