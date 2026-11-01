# 既知三成分上の canonical curve 不変量の保存

2026-11-01。T01、T39=q(5,1)、T40=q(5,2) の一般元について、Poisson 双ベクトル
自体とその倍率を固定し、基点を動かして実際の射影正規化を行う。このとき primitive
一次部分 β から作る canonical Pfaffian curve の **c₄,c₆ はともに一定** である。
j だけの一定性より強い。ここでいう一定性は、固定した一つの双ベクトルの基点に関する
一定性であり、成分全体にわたり同じ値を取るという意味ではない。

この証明は既知三成分に限られる。全S61の未知成分には適用しない。

## 1. 二つの構成に共通する character

五次Pfaffian模型 Φ の変換を

    Φ′(r)=A Φ(Bᵀr) Aᵀ,    χ(A,B)=(det A)² det B

とする。Fisher–Sadek,
[On genus one curves of degree 5 with square-free discriminant](https://arxiv.org/html/1406.6953v1),
§1、§2 Theorem 2.1 の規約では c₄(Φ′)=χ⁴c₄(Φ)、c₆(Φ′)=χ⁶c₆(Φ) である。

一方、Markarian–Polishchuk,
[Compatible Poisson brackets associated with elliptic curves in G(2,5)](https://arxiv.org/html/2310.18759v2),
Introduction の π₅,₂ と Theorem A を用いる。五つの係数 bivector を外積し、
射影接空間へ送る構成は det²V の twist を持ち、滑らかな Grassmannian 線形切断から
q(5,2) を与える。従って、固定した体積形式でこの twist を自明化した双ベクトル ΠΦ は

    ΠΦ′ = χ(A,B) A_* ΠΦ.                       (1)

を満たす。det B は五つの係数の外積により、(det A)² は twist により現れる。
A はここでは同次空間の線形変換であり、右辺の A_* はその射影変換による押出しである。
この式と Fisher の character が同じ向きであることが重要である。

## 2. Graph 模型の原点での値

V=U⊕Ce₀、J₀₁=J₂₃=1 とし、primitive な Z の上三角成分を

    (Z01,Z02,Z03,Z12,Z13,Z23)=(a,b,c,d,e,−a)

で固定する。q(5,2) の復元ノートで得た graph 模型を

    Ψβ(Z)=[ Z                 −J ellβ(Z) ]
           [ (J ellβ(Z))ᵀ      0         ]

とする。点 e₀ での MP 構成の評価は、五つの係数 bivector を
Λ²(V/Ce₀) に射影してからその外積を取ることで計算できる。従って graph の最終列は
消え、固定された五つの primitive Z の基底だけが残る。この五つは独立なので外積は
非零であり、MP 構成の非退化な exterior pairing によって

    ΠΨβ(0)=κJ,    κ≠0                         (2)

を得る。κ は体積形式と基底順に依存する定数であり、β には依存しない。
非零性を滑らかさの仮定や数値標本から推測する必要はない。

保存コードで用いる canonical 行列は

    Φβ=[ star Z       ellβ(Z) ]
        [−ellβ(Z)ᵀ     0     ].

primitive Z では JZJᵀ=−star Z。A₀=diag(J,−1) とおけば
A₀ Ψβ A₀ᵀ=−Φβ である。さらに変数を −I₅ で変換すると Φβ を得る。
この一定変換の character は −1 なので、偶数weightである c₄,c₆ は
Ψβ と Φβ で同じ値を持つ。

## 3. q(5,2) での基点独立性

固定した一般の q(5,2) 双ベクトル P と、その定数核の固定模型 Φ を取る。
`ELLIPTIC_Q52_RECONSTRUCTION_REVIEW_20261101.md` の明示模型は、Fisher の
Hesse 模型 u₅(1,λ) に変数 (r₃/λ,r₄,r₀,r₁/λ,r₂) を代入したものと成分ごとに一致する。
従って一般 λ でその Grassmannian 線形切断は滑らかな楕円曲線であり、MP の定理が
適用できる。核条件と exterior pairing により ΠΦ は P と点ごとに比例する。
両者が非退化Poissonである開集合では、逆二形式の閉性と η↦η∧Ω の同型性から
比例因子が定数であることが従う。よって ΠΦ=αP で、α≠0 は基点に依存しない。
固定した倍率 c≠0 に対して、基点 p を原点に送り、cP を定数項 J、primitive 一次部分
βp にする**全体の射影変換**を A_p とする。perspective による gauge 固定も含める。
核の空間部分への射影が同型である稠密開集合では、基底変更 B_p を一意に選んで

    Ψβp=[A_p,B_p]Φ

と書ける。正規化後の双ベクトルを Pp=(A_p)_*(cP) と置くと、(1) は

    ΠΨβp = χ(A_p,B_p) (α/c) Pp

となる。原点で (2) と Pp(0)=J を用いれば

    χ(A_p,B_p)=κc/α.                            (3)

右辺は p に依存しない。Fisher の変換則と第2節の一定同値変換により
c₄(Φβp)、c₆(Φβp) はともに一定となる。
Sp₄ frame の別の選択はこれらの値を変えない。これは primitive Plücker 空間の
自然な作用と det S=det(Λ²₀S)=1 から同じ変換則で従う。

## 4. q(5,1) とトーリック成分

`ELLIPTIC_GENERIC_DEGREES_20261101.md` の全パラメータ恒等式は、二次Cremona写像 fλ
に対して、正規化済み一次部分の関係

    b₂(λ,fλ(p),−4c)=S₀ · b₁(λ,p,c),    S₀∈Sp₄

を与える。λ,c を固定して p を動かすと右辺に対応する q(5,2) の倍率も −4c に固定
される。Sp₄ 不変性と第3節より、q(5,1) の c₄,c₆ も一定である。
非線形Cremona同値だけを根拠にせず、既に検証した primitive 一次部分の追加恒等式を使う。

T01 の一般の log-canonical 双ベクトルは、開トーラス上で推移的に作用する射影トーラス
によって双ベクトル自体が保たれる。異なる基点の正規化後の一次部分は Sp₄ 同値であり、
同じ不変量を持つ。これで既知三成分についての主張を得る。

## 5. S61 の微分条件としての範囲

基点移動の普遍的座標式を D_hβ=M_hq+a_h(β) と書く。M_h は定数、a_h は β に関して
高々二次である。従って

    dc₄|β(D_hβ)=0,    dc₆|β(D_hβ)=0,    h=1,…,4

は、既知三成分の稠密開部分で成立し、左辺は多項式なのでその閉包にも成立する。
重み deg β=1、deg q=2 に関してそれぞれweight9、13である。
j の分母を払った条件 (2c₄ dc₆−3c₆ dc₄)(D_hβ)=0 のweightは21。

`PRIMITIVE_CURVE_INVARIANT_GRADIENTS_20261101.md` と対応する証明書は β* での
8本の係数階数が8であることを厳密に確認した。
**無条件の全Poissonイデアルまたはその根基への所属は偽**である。後続ノート
`CURVE_INVARIANT_CONSERVATION_COUNTEREXAMPLE_20261101.md` の大域的二次Poisson反例が示す。
さらに `CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_20261101.md` は Δ≠0 上でも所属を否定する。
これらは既知三成分に限定した本ノートの保存則を否定するものではない。
全一般ファイバーをこの37変数部分空間へ制限することはまだ正当化されない。

独立点検：elliptic_jets 担当が (1)–(3)、graph の評価、およびCremonaの倍率規約を確認した。
一次資料本文の該当箇所は本ノート執筆時に再確認した。
