# nil31長translation第二全五源のT40帰属：独立レビュー

2026-11-07。主義務 `E_reduced_nonnormal`。
**H+既約quarticのnil31長translation型の第二完全五源と、その全PGL飽和はT40に入る。**
全パラメータ域a≠0での実affine同値と、canonical点のfaithful S61全graphにより確認した。
第一源の新成分性・相違は別担当の証明であり、本レビューで重複して検証していない。
本結果で追加する成分数は零である。

独立実装は `code/classification_nil31_long_second_review.py`、
報告は `reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_SECOND_REVIEW_20261107.json`、
有限契約は `jobs/hyperplane_quartic_nil31_long_second_review.json`。
生産者、そのhelper、S61算術engine、CAS、NumPy、legacy算術をimport・実行しない。
既存の別作者Fraction/Laurent算術と、旧独立source reviewerの基本関数だけを再利用した。
最終実行5.716秒、最大RSS70,924KiBで全検査PASS。
既済producer・証明書のbytesと、freezeされた100頁主TeXを変更していない。

## 1. 元90係数から第二完全五源を保つ

全nilpotent源のrow4について、divergenceとLieの全110線形行を別に再構成した。
その完全解は15係数であり、表示源が全行を満たして係数写像rank15を持つことを確認した。
したがって選んだ部分族から計算を始めていない。

第二全源の自由係数は `(a,n3,n4,n6,n12)∈C*×A4`、a=n14である。
独立に入力した全15置換は

    n14=a, n13=−1/6, n2=1/(216a²), n8=−n12/(6a),
    n1=−27a n4²−n6/(12a)−n12/(36a²),
    n0=54a n4²n12−18a n3n4+n6n12/(2a),
    n5=n12²/a−6n6n12, n9=1/(36a),
    n11=−9a n4, n7=−6a n6+2n12,
    n10=6a(3n4n12−n3),
    n3,n4,n6,n12 は自由。

これは元の支持証明書の全15式と一致し、元の六括弧へ戻した全係数も一致する。
全五変数を残したまま十個のstrict homogeneous Jacobiが恒等零である。
a以外の源パラメータは一度も反転しない。
全a=0補集合の排除と「第二五源が全残支持の一つである」という完全性は、
別に独立受理済みの `HYPERPLANE_QUARTIC_NIL31_LONG_SUPPORT.md` と
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_REVIEW_20261107.json` を依存として用いる。

## 2. 全五源を一つの実PGL軌道へ送る

canonical点はa=1,n3=n4=n6=n12=0、すなわち

    n2=1/216,n9=1/36,n13=−1/6,n14=1、他の元15係数は零。

d=−3n6/2+n12/(2a)、s=18a n6、f=−3a n4、e=−2n3−54a n4n6+6n4n12と置く。
順序 `(h,x,y,z,w)` の行列Mは

    h'=h,
    x'=x+d z+e h+f w,
    y'=y+d h,
    z'=z,
    w'=a w+s h。

その行列式はaであり、逆は

    h=h', z=z', y=y'−d h',
    w=(w'−s h')/a,
    x=x'−d z'−(f/a)w'+(fs/a−e)h'。

全行列係数をQQ[a,a⁻¹,n3,n4,n6,n12]で独立に作り、MM⁻¹=M⁻¹M=Iを全成分確認した。
またhomogeneous cofactor `N=y∂x+z∂y+h∂z` との可換性も全行列で成立する。
canonical tensorをP0、普遍五源をPとすると、全十個の添字対の全多項式係数で

    P(MX)=M P0(X) Mᵀ

を検査した。tensorのscalar係数は1である。
したがってa非零の全域、自由係数が零となる全境界を含め、各点はP0の実affine/PGL像である。
稠密軌道の推測、Pfaffianの因子型だけの一致、formal jetの一致を使わない。

## 3. canonical点の全周囲接次元25

P0を216倍した整数tensorから、全150変数の全350射影Jacobi微分Dを独立に作った。
全有理消去と独立minorでrankD=100。
producerに保存されたliteral100次minorも元の行・列で別に再計算して一致した。
保存された全50本のkernelを全150係数へ復号し、全350行でDK=0、rankK=50を確認した。
よってrankの下界と上界が一致し、選択部分行列だけの接空間ではない。

Eulerの全25列と、全GL5列、tensor scalar列も別に作り、D積が零であることを確認した。
Euler rank24、群・垂直・scalarを合わせたrank47である。
従って射影全Poisson空間の接次元は

    (150−100)−24−1=25。

実PGL軌道は22次元なので、このorbitだけから成分性を推論することはできない。
以下で既知25次元T40を通すことが必要な橋になる。

## 4. 全150商の正逆、全134式と全45T40 graph

P0のaffine評価 `(x,y,z,w)=(0,1,0,0)` の実Pfaffianは正確に4/3。
y'=y−hの実translationを独立に施した原点での全反対称定数行列は

    [ 0   −1/3   0    2 ]
    [ 1/3   0    1/3  0 ]
    [ 0   −1/3   0   −2 ]
    [−2     0    2    0 ]。

保存された有理Darboux frame BについてBKBᵀ=Jを全成分確認した。
続くperspectiveは `(0,−5/18,0,1/6)` であり、translation、B、perspectiveの積Mと
その逆を有理行列で独立構成し、両積が単位であることを確認した。
反転するのは非零の有理定数だけで、源parameterの新たな例外はない。

全61個のS61値を元のraw free indexに入れ、旧modelの全消去変数式から全150係数を再構成した。
実際に変換したtensorとの差をΔとすると、単に射影外積を零と確認するだけでなく

    Δ=Euler∧div(Δ)/5

を全150係数で確認し、明示されたlinear Euler gaugeを保存した。
再構成したtensorを逆PGLで戻した側でも、同じ全150係数のEuler gauge等式を確認した。
したがって特殊点を部分的な座標へ投影しただけではなく、faithfulな正逆商同値である。
元の全134方程式も全て零となった。

既存 `Q52_UNIVERSAL_SECTION_20261101.json` の全45 graph式を、全61値で独立評価すると全て零。
Q51/T39の全45式は零とはならず、その全残差も保存されたものと一致する。
どちらも一部のhead/graph式だけを検査したのではない。

ここで使う既済外部入力は「普遍Q52 graphは既知25次元T40成分の被約S61閉包そのもの」という
全graph認識定理である。この一つの点の代入でその普遍定理を新たに証明したとはしない。
正逆faithful商、全134式、全graphにより、P0は実際にT40に入る。

## 5. 正則唯一性と全源への結論

既知T40の次元25がcanonical点の局所次元の下界を与え、第3節の埋込み接次元25が上界を与える。
従って局所環は正則、同点を通る全空間成分は唯一で、それはT40である。
さらに第2節の全a非零域の実同値と、連結PGLが既約成分を保つことから、
第二完全五源の全点と全射影飽和がT40へ入る。
源の一般点だけを閉包に入れて零境界を推測する論法は不要である。

本レビューは第二源を全被覆し、新成分を追加しない。
型5全体のE完了には、別担当の第一源の全周囲成分性・既存成分との相違、
及び既済の全境界排除を合わせる必要がある。
G全体、既約非正規五次全体、主分類完成について新たな結論を加えていない。
