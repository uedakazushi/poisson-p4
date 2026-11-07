# dependent globalV全Jordan源は五fixed＋一parameter線へ縮む

2026-11-07、主義務G。producer手証明。
normalss部分はrootが別に導出した係数式と照合した。
normalnil trace零部分は既存別作者review済みnoteに一致する。
全higher support又は残六源のG非優越はまだ本noteから結論しない。

## 1. 二つのnormal hyperplanesがPoissonになる補集合

被約既約actual quintic Dとrank零planeを仮定し、新modular入口のfull35/30/9 sourceを使う。
traceless common field V(x)=A*x≠0、trA=0、trace field T=tVを保つ。
full Π1には`w² C∧∂z`と`z² Dv∧∂w`があり、C,Dvはconstant-base vectors。

normalssではK=kzw、normalnilではK=αz²であり、z=0でnormal contractionに残るのは符号を除きw²Cだけ。
normalssのw=0でも同様にz²Dvだけ。
従ってC=0、又はnormalssでDv=0なら、そのlinear hyperplaneがPoissonとなりactual退化因子にlinear factorが現れる。
GのDは既約quinticなので、normalnilでC≠0、normalssではC,Dvの両方非零が必要。
Π2は全base-baseなのでこの判定を変えない。

## 2. 三次元vector spaceの一様なwedge lemma

C≠0とし、全xについて

`x∧(pAC+kC)−q(Ax)∧C=0`、`q≠0`。                      (1)

p≠0ならquotient W/Cへ射影すると、全xbarで`xbar∧(AC)bar=0`。
quotientは二次元なのでAC∈C-line。
p=0でも(1)はim(kI−qA)⊂C-lineを直接与える。
ゆえに両場合ともあるcとlinear form ℓで

`A=cI+C⊗ℓ`、`ℓ(C)=−3c`、`AC=−2cC`。

(1)に戻すと`k=(2p+q)c`である。
c非零ならAはdouble-eigenvalue semisimple型、c零ならA≠0からrank一nilpotent型。
この証明はx=Cから得る補助係数の例外に依存せず、q非零の全値を覆う。

## 3. normalss：全Jordan型からH固定とN線だけが残る

全normalss symbolは

`Π0=V∧((t+1)z∂z+(t−1)w∂w)+kzw∂z∧∂w`。

`a=t+1`, `b=t−1`と置くと、full projective Jacobiのw³ coefficientは

`x∧(bAC+kC)−(2b−a)(Ax)∧C`。                         (2)

これはstrict成分`Jijz=(2b−a)w²(V∧C)ij`、`Jijw=0`、`Jizw=w³(bAC+kC)i`から得られる。
他の全Π1項はnormal monomialが異なり、full higher termsはweight二以上なのでこの係数へ寄与しない。
z³ coefficientも全sourceで

`x∧(aA Dv−kDv)−(2a−b)(Ax)∧Dv`。                      (3)

normal permutation z↔wで(t,V,k)→(−t,−V,−k)とするactual逆変換からも(3)を得る。

(2)のq=t−3、(3)のq=t+3は同時に零にはならない。
C,Dv非零と§2から、全tでVはdouble semisimple又はrank一nilだけ。
regular semisimple、mixed、regular nilpotentの全sourceはG空。
三parameter sourceをgeneric periodだけで除外したのではなく、全exception t=±3も一方の式で覆う。

V=H=diag(1,1,−2)に直す。
t≠3ではCはdistinguished eigenlineに入り、(2)はk=3t−5。
t=3では(2)は(2A+kI)C=0となりk=4又は−2だけ。
同様にt≠−3ではDvがdistinguished line、k=−3t−5。
t=−3の例外はk=4又は−2。
両式を満たすには**t=0,k=−5だけ**となる。
t=3で(3)が要求するk=−14、t=−3で(2)が要求するk=−14は例外の4/−2と一致しない。

V=N=E01に直すと、q非零の側でc=0なのでk=0が必要。
これによりnormalssの全sourceで残すwhole族は

- `V=H,T=0,K=−5zw`、C,Dvはdistinguished lineの非零vectors。
- `V=N,T=tN,K=0`、t全C。

である。後者のt=±3で許されるC又はDvの拡大も保持する。

## 4. normalnil：trace零の二残余とtrace非零のN二orbit

全normalnil dependent symbolは

`Π0=V∧(t Euler_normal+z∂w)+αz²∂z∧∂w`。

t=0の全証明は`PLANE_RANK_ZERO_DEPENDENT_TRACE_ZERO_HYPERPLANE_REDUCTION.md`にあり、
別作者reviewを受けた。
zw² coefficientは`x∧(AC+2αC)−3(Ax)∧C`。
§2でもc=2α/5となるので、残る全sourceは

- `V=H,T=0,K=(5/2)z²`。
- `V=N,T=0,K=0`。

だけである。

t非零ではfull weight一Jacobiのw³ coefficientが

`t*(x∧AC−(Ax)∧C)`。                                (4)

normalnil非対角部分及びαz²の寄与にはz因子があるので、(4)へ混ざらない。
他のΠ1項も同じ極端normal monomialへ寄与しない。
§2にp=q=t,k=0を適用するとc=0、従ってVはrank一nilpotent型だけ。
actual independent trace/common-field scalarsでT=V=Nへ直せる。
base nilpotent gradingとtensor scalarでα零/非零は0/1へ直り、両whole orbitを保持する。

**ここからα零とは結論しない。** t非零のzw² coefficientには他のΠ1項も寄与するため、
t零の係数式をそのまま適用してα非零を落とすことは出来ない。

## 5. 残六whole源と有限target

全actual normalizations及び逆はdependent first-jet noteと上のscalar/base/normal changesにより得られる。
全二十一familyからGで未処理として残すのは次の六だけである。

| model | trace T | common V | K | whole parameters | affine symbol base B |
|---|---|---|---|---|---:|
| ss H | 0 | H | −5zw | fixed whole source | 7 |
| ss N | tN | N | 0 | t全C | 7、fixed tでは6 |
| nil H | 0 | H | (5/2)z² | fixed whole source | 6 |
| nil N trace零 | 0 | N | 0 | fixed whole source | 5 |
| nil N trace非零 | N | N | z² | fixed whole source | 7 |
| nil N trace非零/K零 | N | N | 0 | fixed whole source | 6 |

Bはdependent first-jet noteの全source countsを戻した値で、選んだanchorのdimensionではない。
full higher kernel sumをdとすれば、plane六とprojectivization一を戻すimage upperはB+d+5。
必要ならGのC≠0、normalssではDv≠0も全necessary source条件として保持出来る。
まず全30/9と全t exceptionを一つのfinite jobで処理し、nonclosingならfull original Jacobi又は実Pf因子構造へ進む。
この還元だけで残六sourceのG非優越又は全高m完了を宣言しない。
