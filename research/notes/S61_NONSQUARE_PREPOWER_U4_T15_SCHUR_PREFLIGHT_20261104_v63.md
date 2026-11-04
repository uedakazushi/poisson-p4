# direct P4+H5 の一つの実 Schur 成分：源・support・高さ・lifetime

対象は compact H 源u=4（current15行5）、元32標的15、Mark の direct residual6×13 の位置(0,9)だけである。全行や他の成分を計算しない。

受理済み P4 は源 [Q9,G18,G39,G68]、標的 [22,11,13,25]、行列式 D4=a^4 C^3 R434。H5の先頭pairは compact 行[1,3]・列[0,4]、行列式 delta2=-(6/5)B p²。u=4の残る H 三成分は literal zero なので、この一源を H5 で消すのに lower3 の逆行列は必要ない。

rhoL を rho_u,rho_pair の既知単元 monomial LCM（定数因子は1とする）とし

    Bhat=rhoL*(delta2 P_u/rho_u
                -(N_u,pair Adj(A2))*(P_pair/rho_pair))

と置く。欲しい cleared 成分は

    Saux=D4 Bhat[15]-Bhat[piv4] Adj(P4) P4[15]。

Mark の HYBRID Schur 成分は厳密に

    rho_u/(delta2*rhoL*D4) * Saux

である。さらに Saux の既知 a Laurent 冪とB分母冪だけを払えば多項式 numerator となる。h22 の逆元は使わず、この一成分に Delta3 の逆元も現れない。全六行の chart は引き続き FULLU∩D(Delta5 R434) であり、p/Rやその補集合を捨てない。

## 実値を計算しない事前監査

前の固定3小行列式に保存された全108値を可能な限り再利用する。追加の scalar は既存 Q2L DAG の source atom だけである。P4の全十六実入力と既存818項D4を使って、随伴行列は3×3 cofactorを一つずつ構成・消費する。十六cofactorを一度に保持する方式にしない。

型付き演算DAGに、Laurent a指数箱・B分母冪・support上界、共通QQ分母Dと整数係数絶対値和上界Nを付ける。和ではB分母の最大冪に揃え、Bの係数絶対値和179を使う。積・和の整数majorantだけを計算し、パラメータ多項式係数は掛けない。共通QQ分母の上界を勝手に約分しない。

最後の一成分から必要なDAG祖先を取り、各値の最終使用回数によりlive項数を数える。演算中の追加output-size scratch、literal入力、15000項の一般作業reserveも加える。support34650・live126957という初期BOX見積りをそのまま仮定せず、これらの再現可能な細分化した上界で実行可否を判定する。32MiB/150000項ゲートが不合格なら実成分計算には進まない。

対応する `operator-and-cost.json` は全source参照、演算、分母、clearing係数、Markの不変な演算根との対応を保存する。実成分は零や非単元になり得る。非零一成分でも全六行の階数や global full57 は証明しない。
