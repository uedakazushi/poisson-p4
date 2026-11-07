# 滑らかなprojected Veronese rank零曲面枝の独立全排除レビュー

2026-11-07。主義務 `E_reduced_nonnormal`。
`PROJECTED_VERONESE_ZERO_SURFACE_EXCLUSION.md` の結論を受理する。
**複素P4の滑らかな次数四projected Veronese曲面上で零となるPoisson sectionは零のみである。**
rank四だけでなく、非零のrank二Poisson sectionも存在しない。
他の次数四曲面、特異射影、他のrank零scheme、既約非正規五次全体を分類したという意味ではない。
正規化Xの滑らかさ、ordinary-double型、導手の滑らかさは本証明の仮定にない。

独立実装 `code/classification_projected_veronese_zero_review.py` は標準ライブラリの
`Fraction` と独自の疎多項式・有理行消去だけを用いる。
producer、そのhelper、CAS、NumPy、legacy、既存リポジトリ算術を一切importしない。
producer JSONは係数比較用のデータとしてのみ読む。
報告 `reports/CLASSIFICATION_PROJECTED_VERONESE_ZERO_REVIEW_20261107.json` は
実測0.297秒、最大RSS17368KiB、全検査PASSである。
有限契約は `jobs/projected_veronese_zero_review.json`。gateは数学証明の代用ではない。

## 1. 全軌道を一つの埋込みへ還元する幾何

標準Veroneseの座標を `(a²,b²,c²,ab,ac,bc)` とすると、線形写像

    (u0,u1,u2,u3,u4,u5) ↦ (u0−u1,u1−u2,u3,u4,u5)

のkernelは `(1,1,1,0,0,0)` の直線である。
対応する対称行列は単位行列、従ってrank三。
Veroneseのsecant varietyはrank二以下の対称行列全体で、tangent方向もその中にある。
この中心からの射影は点と接ベクトルを分離し、properな閉埋込みとなる。
その像が

    S=[a²−b²:b²−c²:ab:ac:bc]

である。全rank三対称行列は複素数上GL3の合同変換で同値であり、
Sym²表現はVeroneseと射影中心を同時に移す。
従って商五次元空間上のPGL5変換を誘導し、全非secant中心の射影は一軌道である。

次数四の滑らかなprojected VeroneseをP2のO(2)の五次元部分系として表すと、
その中心はrank三でなければならない。
rank一中心はVeronese自身の点であり、この次数四の閉埋込みを与えない。
rank二中心は相異なる二点を結ぶsecant上にある。
このとき射影はO(2)がampleなので有限であり、中心のrank二supportに対応するP2内の
固定直線の外では一対一であるため双有理であるが、secantの二点を同一視する。
像が滑らかなら有限双有理射は同型となるため矛盾する。
従って上の一軌道で、今回の「滑らかな次数四projected Veronese」の全てを尽くす。
一般の次数四曲面をこの軌道へ帰着したわけではない。

## 2. 完全な線形源とEuler商を独立構成

全150係数ΠのsectionがS上で零であることは、各点のEuler商で
`E∧Π=0` と同値である。Sの五個の二次式をq_iとすると、全十個の三添字について

    q_i Π_jk(q) − q_j Π_ik(q) + q_k Π_ij(q)=0

を課す。各式の全28個の六次単項式係数を展開し、全280×150行列Aをゼロから再構成した。
保存された全行・全列と正確に一致し、独立有理RREFから

    rank A=121, dim ker A=29

を得た。literalな121行・121列のminorも独立計算し、行列式は正確に1。
29本の独立kernelを独立に生成し、全行でA・kernel=0を確認した。
producerの29本も全てAで消え、独立であることを別に確認した。

Euler∧linear場の全25列を別に構成し、全Aで消えること、階数24を確認した。
全29本を `Π−E∧div(Π)/5` でcanonical化し、その行空間を独立RREFすると五基底を得る。
五基底はdivergence零で、保存された全十括弧と一致する。
Euler列とこの五列の和の階数は29であり、従って選んだ部分族ではなく完全商である。
Euler列の外積完全列とH1(T_P4)=0により、これが
`H0(Λ²T_P4⊗I_S)` 全体を与える。
五係数z0,…,z4を自由なまま保持する。係数の逆元や非零chart条件はない。

## 3. 全350Jacobiと根基原点の明示証明

全五源の実括弧からSchouten三ベクトルを独立展開し、Euler外積を取り、
五個の四添字と全70個の四次単項式の全350方程式を得た。
保存された全350式と正確に一致する。
二次単項式十五個に対する全350×15行列の階数は14であり、次の十四式が同じ行空間を張る。

    r0=z0²−z4², r1=z0z1−z4², r2=z1²−4z4²,
    r3=z0z2, r4=z1z2, r5=z2²−z4²,
    r6=z0z3, r7=z1z3, r8=z2z3, r9=z3²−z4²,
    r10=z0z4, r11=z1z4, r12=z2z4, r13=z3z4.

各r_iを元の全350行で表す有理functionalを、独立行消去から生成してliteralに検査した。
逆に各350式を十四式で表す係数も独立に生成・検査した。
両向きの全係数をレビュー報告に保存してあり、全Jacobi idealが `(r0,…,r13)` である。

次の五つの恒等式を疎多項式として独立に検査した。

    z0³=z0 r0+z4 r10,
    z1³=z1 r2+4z4 r11,
    z2³=z2 r5+z4 r12,
    z3³=z3 r9+z4 r13,
    z4³=z0 r10−z4 r0.

各変数の三乗がJacobi idealに入り、全式は原点で零だから、根基は正確に
`(z0,z1,z2,z3,z4)`。従って被約Poisson支持は原点だけである。
有限jet、標本点、数値rank、Groebner black box、flatnessを使っていない。
Euler方向だけのliftは射影sectionとして零なので、別の非零Poissonを与えない。

## 4. Pfaffian規約と結論の範囲

全五個の符号付き4×4 Pfaffian minorも独立に作り、producerの全係数と一致を確認した。
正しい射影anticanonical五次式は

    D=Σ_i (−1)^i x_i Pf(Π_without_i)

であり、これは全五源上でJacobiを課す前から恒等零である。
別にx0=1の実Euler補正代表のPfaffianも独立に計算し、恒等零を確認した。
従ってrank四不在だけなら線形源全体の時点で既に従う。
非零Poisson自体の不在には前節の全Jacobi根基を使う。

初回producerの「五個のminorがx_iと共通五次式の積になるか」を探索する診断は、
一般のhomogeneous liftについて正しい恒等式ではない。
それを一切使わず、上のEuler収縮と実affine代表で正しい式を確認した。
producer旧bytesを修正したのではなく、明示的にその診断を証明から除いている。

PGLはPoisson条件とsectionの零を保つので、第一節の一軌道性から全滑らかなprojected
Veronese曲面上の非零Poisson sectionを排除する。
これだけで他の非ACM零曲面や全既約非正規五次を排除したとはしない。
新成分数や全分類の完成度を変更する結果でもない。
