# 全射影NのSerre R1：階数二境界因子の滑らかな一点

2026年11月3日。これは全射影27次元の被約Nについて、S61の外も含めR1を証明する追補である。
全正規性・Cohen–Macaulay性は主張しない。指定Q_(1)での非正規・非CM性はそのままである。

Nのproper模型は、9次元旗上の不変19次元空間f∈hWのP18束である。
全rank≤2境界はNg(dim26)とNf(dim25)の二成分。
rank4開のR1は既に証明したので、唯一のrank2境界因子Ngの一点でNが滑らかと示せば十分である。

## 明示点と全源の微分

基準旗L=C∂y、W=span(h,x)、H=Chを取り、

    f=h*x、g=0、potential=z³+w³

とする。生の構造はP=∂y∧Vで、

    V=−hx∂x−3w²∂z+3z²∂w

である。canonical代表は∂y∧C、5Cは

    Ctilde=(h²,−4hx,−15w²+hz,15z²+hw)

となる。

proper模型の全19係数方向とGL5作用をcanonical150係数へ送った全affine微分の階数は28。
非零小行列式6912/390625を含む全行列を保存した。
源はprojective27次元なので、その全微分はこの点で単射である。
ここでまだtargetの接空間と源の微分像を同一視しない。

## proper模型の全幾何学的ファイバーは一点

Eul∧v∧P_can=0の全175係数行列は階数4で、∂yが核に入るので、定ベクトル因子は一意。
Ctildeの四成分は独立な二次式なので、P_canは定covector Casimirを持たない。
一方Nfの全点は二つの独立な定Casimirを持つため、この点はNfに属さない。
したがって任意のproper模型の持上げはNg側で、同じ定ベクトル因子Lを持つ。

Ngでの旗には、商C5/L上の独立な一次式ell,kと定数lambda,alphaが必要で、

    Ctilde(ell)=lambda*ell²,
    Ctilde(k)=−alpha*ell²−4lambda*ell*k

を満たす。これはcanonical化によりEulerの曖昧さを除いた内在条件である。
ell=A h+B x+C z+D wとする。
lambda=0なら四成分の独立性からell=0。
lambda≠0なら、第一式のx²係数は−lambda B²でB=0。
zw係数からCD=0、z²,w²係数から

    15D=lambda C²、−15C=lambda D²

となりC=D=0を得る。従ってellはhの倍数で、lambda*ell=hとなる。

第二式はCtilde(k)+4h kがh²の倍数という条件である。
k=A h+B x+C z+D wに対して左辺は正確に

    5A h²+C(−15w²+5hz)+D(15z²+5hw)

なのでC=D=0。kはhと独立だからB≠0で、W=span(h,x)に一意に決まる。
こうしてL,W,Hという旗全体が戻る。
固定旗ではcanonical20係数写像が単射なので、テンソルの係数も一意。
proper模型の**全幾何ファイバー**は一点である。

## 長さ一からtargetの滑らかさへ

全源の微分が単射なので、上の一点ファイバーの接空間は零である。
幾何的台が一点、有限型で接空間が零だから、schemeファイバーは被約な長さ一。
propernessによりこの点の近傍で射は有限となる。

target環Aから有限な源環Bへの写像は閉像Nを定める単射であり、
B/m_A B=Cを満たす。有限加群B/AにNakayamaを適用するとB=Aとなる。
従って射はこの点の近傍で同型で、滑らかな源からN自身が滑らかと分かる。
単なる源の微分の階数や、連結ファイバーだけからtargetを滑らかとしたのではない。

Ngは既約な26次元因子なので、そのgeneric pointは滑らかである。
Nfは余次元2、rank4開ではR1なので、全Nの全余次元1点は正則。
これで全射影NのSerre R1が証明された。

## 証拠

- `code/ambient_N_global_R1_20261103_v63.py`
- `reports/AMBIENT_N_GLOBAL_R1_20261103_v63/report.json`

全QQ行列、minor、内在flag条件の係数恒等式を約0.7秒で検証した。
既存のglobal proper N模型と、26/25次元の全境界被覆を明示的に使う。
全正規化、全境界の特異点分類、全交叉scheme構造は別課題である。
