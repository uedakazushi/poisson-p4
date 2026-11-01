# AB の原点：三次の正規化関数も一次座標で生成される

2026年11月1日。B=C[AB]、Dをその正規化とする。
既証明のproper birational incidence X→ABの正規性からD=H⁰(X,O_X)である。
XのLG(2,4)への射影に沿う次数n部分をE_nと書けばD_n=H⁰(LG,E_n)となる。
[二次の結果](AB_ORIGIN_DEGREE_TWO_20261101.md)を三次へ進め、

    dim D₃=dim B₃=6691,  H¹(LG,E₃)=0

を証明する。従って次数三までに追加の正規化生成元はない。
これは全次数の生成性や原点の正規性の証明ではない。

## 1. 四段階filtrationのcohomology

固定L=span(∂y,∂w)、q=(x,z)、p=(y,w)を用いる。
固定平面上のincidence fiberはY×A³で、Yはsl₃内の階数一以下の円錐である。
その座標は(s0,s1,s2,h0,…,h4)、垂直二次式はg0 x²+g1 xz+g2 z²とする。
正規化されたN(s,h)とQへの線形写像は、二次ノートと同じである。

垂直座標gの次数によるP-stable filtrationを取る。unipotent radicalはgをNの線形式だけずらすので、
このfiltrationは分裂直和とは限らない。三次のassociated gradedは

    Y₃ ⊕ (Y₂⊗V₁) ⊕ (Y₁⊗V₂) ⊕ V₃

である。Y_nの次元は(n+1)³。2×2 minorsのGroebner基底の標準単項式から各Levi characterを計算し、
Borel–Weil–Bottを適用すると、associated gradedのcohomologyはH⁰,H¹に限られ、

    dim H⁰(gr E₃)=7285,
    H¹(gr E₃)=W_(9,1) ⊕ W_(8,0),
    dim W_(9,1)=429,  dim W_(8,0)=165

となる。W_(a,b)はSp4のpartition表記であり、W_(8,0)=Sym⁸Vである。
H⁰(gr E₃)内のこれらの重複数は、それぞれ一、二である。
全Levi分解とBott表を証明書に保存した。

Bottの規約は二次ノートと同じで、U*のweightにρ=(2,1)を加え、符号付き置換でdominant chamberへ移す。
一次資料は[Fonarev, Theorem 3.2](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9880&wshow=paper)であり、
同論文のUからU*への双対weight変換は既存の独立reviewで点検済みである。
ここではcharacterのEuler差だけを最初に計算しており、594次元のH¹が消えることは次節で別に証明する。

## 2. Parabolic評価からH¹の消滅を強制する

PはLを保つparabolicである。点のpushforwardによるfiber関数のpullback規約を全ての表現で統一する。
Leviの対角行列に対し

    wt(h0)=(-3,1), wt(h1)=wt(s0)=(-2,0),
    wt(g0)=(-3,-1), wt(g1)=(-2,-2)

である。Pは、pに含まれるlowest vectorからq方向へ上げる正根部分群とLeviを含む。
従って各既約Sp4表現のlowest vectorはP-cyclicである。
具体的にはSym⁸Vではy⁸、W_(9,1)ではCartan成分内の

    y⁸ (y∧w) ∈ Sym⁸V ⊗ Λ²₀V

を取れる。後者のweightは(-9,-1)である。
このcyclicityは、Pがlowest vectorに対する反対Borelを含み、その正根部分群が
既約表現をlowest vectorから生成することから従う。

非零のG-submoduleをH⁰(E₃)に取ると、固定平面での評価はそのsubmodule上で零にはならない。
もし零ならGの平面に対する推移性により全ての評価が零となるからである。
従ってFrobenius reciprocity、またはこの評価写像により、各G型の重複数は対応するP-map空間の次元で抑えられる。
P-cyclicityによりP-mapはlowest vectorの像で一意に決まる。

### W_(9,1) は現れない

E₃内のweight(-9,-1)はh0 g0²が張る一次元空間である。
radicalのshear w→w+u zは、p→p+Cqの対称行列CでC[1,1]=u（添字は零始まり）、他を零としたものである。このshearはNを固定し、

    g0 → g0+12u h0

と作用する。一方、y⁸(y∧w)にはuの二次項がない。
よってP-mapの像c h0 g0²には二階微分零が必要だが、

    (d/du)²(c h0(g0+12u h0)²)=288c h0³

である。h0³はY上で非零であり、例えばN=E13で非零となる。
従ってc=0で、cyclicityからP-map全体が零である。
H⁰(E₃)にW_(9,1)は含まれない。

### Sym⁸V の重複数は高々一

E₃内のweight(-8,0)は次の三つの独立な関数が張る。

    h0²g1,  h0h1g0,  h0s0g0.

Yの稠密なchart N=(1,a,b)ᵀ(-ad-be,d,e)で独立性を直接確認した。
lowest vector y⁸の像を

    f=A h0²g1+B h0h1g0+C h0s0g0

と置く。Levi upper shearはyを固定し、h0,s0,g0を固定する一方、

    h1 → h1−4u h0,  g1 → g1−2u g0

と作用する。従って2A+4B=0が必要である。
またw→w+u zはyを固定し、

    g1 → g1+u(6h1−s0/2)

と作用するから、6A+12B=0と−A/2+12C=0が必要となる。
これらの条件はY上でも階数二であり、可能な像は

    h0²g1 − (1/2)h0h1g0 + (1/24)h0s0g0

の高々一次元空間に限られる。P-mapの存在をこれらの必要条件だけから推論する必要はない。
重複数の上界一が得られれば十分である。

## 3. Euler characterと正規化次元

有限filtrationのassociated gradedにはH²以上がないので、E₃にもH²以上はない。
H¹(E₃)に現れ得るG型はW_(9,1)、W_(8,0)だけである。
表現環内のEuler characterはfiltrationで不変である。

W_(9,1)のEuler重複数は1−1=0で、H⁰重複数の上界は零であるから、H⁰,H¹とも零。
W_(8,0)のEuler重複数は2−1=1で、H⁰重複数の上界は一であるから、H⁰は一、H¹は零となる。
他の型にはH¹がない。従ってH¹(E₃)=0、

    dim D₃=7285−429−165=6691

を得る。上のSym⁸候補を見ただけで165次元の追加生成元が存在すると結論してはいけない。
associated gradedにはその型が二コピーあり、一コピーが残ることはEuler差と両立する。

## 4. 座標三次式の厳密な下界

LGのbig cell三変数、Yの稠密chart四変数、垂直g三変数から成る十変数のincidence chartを使う。
そこへQの60個の係数関数をpullbackし、45個の独立座標を取る。
その全16215個の三次積を疎多項式として展開し、素数32003で消去した。
6691個の積と6691個の単項式で定まる小行列のdeterminantは、三角化pivot積として

    20093 mod 32003 ≠ 0

である。全分母は32003で可逆である。従ってQQ上の階数は6691以上となる。
B₃⊂D₃と前節の幾何学的上界6691から等号が従う。
有限体での階数を標数零の上界として用いてはいない。

## 再現と未決事項

```
.venv/bin/python code/ab_origin_degree_three_20261101.py --output /tmp/ab-origin-degree-three.json
```

証明書はBott表、parabolic作用の必要条件、6691次小行列の添字とpivot、入力hashを保存する。
B₀=D₀、B₁=D₁、B₂=D₂と合わせ、次数三まで正規化関数は一次座標から生成される。
原点の全次数生成性、正規性、導手idealと実際のT01零ファイバーのスキーム構造は未決のままである。
