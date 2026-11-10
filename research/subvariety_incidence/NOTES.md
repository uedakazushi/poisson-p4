# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v7。研究方針は RESEARCH_POLICY.md、全証明は現行 papers/poisson_subvariety_incidence_ja.tex と共有TeXにある。

## 共鳴Logの新しい等方的境界

$V=U\oplus W$（dim U=3, dim W=2）、重み(1,1,-2)のLog源を$\gamma=\iota_\xi(du_0\wedge du_1\wedge dz)$、$\iota_{P_0}\Omega=dq\wedge\gamma$で作る。従来の直交分解$q=q_U+w_0^2+w_1^2$の$q_U$退化は$q$自体の退化である。滑らかなQを保つ境界は$q|_W$の退化で調べる必要がある。

三つの確定した全延長ファイバー：
- q1=u0²+u1²+z²+2u0w0+w1²、rank(q|W)=1：全延長A²_(a,b)、全Jacobi ideal 0、rank4ならa≠0、Pfaffian実因子Q+2H(u0)+H(z)。
- q2=u0²+2u1w0+2zw1+w0²、rank(q|W)=1：全延長Spec C[a,b]/(a(b−4))、rank4枝b=4,a≠0、因子Q+H(z)+Q2((u1+w0)²+2zw1)。
- q3=u0²+2u1w0+2zw1、rank(q|W)=0：全延長Spec C[a,b,c]/(ac,c²)、被約台A²_(a,b)、a≠0でrank4、因子Q+2H(u1)+H(z)。nilradical (c)はa=c=0に沿う非零の平方零方向を持つ。被約147成分の追加を意味しない。

qの半Hessianの行列式は(-1,1,1)、strict stabilizerの係数行列rankは(8,8,7)、minorは(2304,-2304,3456)。これを発散Jacobiの必要条件と合わせ全延長候補に還元し、正確な全Jacobiイデアルを計算する。生成コードとJSONはそれぞれcode/subvariety_incidence/check_isotropic_log_boundaries.py、research/subvariety_incidence/certificates/isotropic_log_boundaries.json。既存の8 producer/certificateと主稿を同じGitHub mainから再現する。

## 現在残る数学的課題

rank(q|W)単独ではファイバーを特徴づけられない。rank1の二例がA²とノードを与えるため、必要なのはLog重み空間と二次形式の混合配置・Lie stabilizerとJacobi障害である。特に非被約の例に現れるembedded方向の一次・二次障害を、既存のambient147成分のscheme局所構造と比較することを優先する。

次点：special Aff、4次del Pezzo例外層、P3横断対称性、147成分の幾何学的特徴づけ。新規性の一次資料による包括確認は未完了。既存一般Q3三型、smooth cubic、del Pezzo一般剛性の成果は主稿から再取得できる。

GitHub main正本、TeXと正確な入力・証明書のみ追跡、PDF非追跡、Actionsなし、Libraryは実push失敗時の緊急退避のみ。
