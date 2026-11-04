# B minor22本からの低コスト互いに素なpair選択

凍結済み22本のみをdについて整理し、Sylvester行列次数、a/e次数上限、項数の
順に定めたpairをexact native GCDで調べる。最初の互いに素なpairで停止する。
全試行のliteral quotientと再構成、選択pairの別実装GCDも保存する。
終結式、GB、因数分解、点や有限体での評価は計算しない。

選択pairのQQ[a,e,d] GCDが1ならGauss補題によりQQ(a,e)[d]でも互いに素であり、
d終結式は非零である。実際のB-rankbad集合はそのpairの共通零点に含まれ、
そのa/e射影は終結式の零点に含まれる。ただしこのpreflightは終結式そのものも
有限射影も証明しない。有限freeモデルを作るにはdの最高係数消失の場合を別に
覆う必要があり、genericなd次数を全域で固定したとみなさない。

係数上限は各QQ多項式を分母払後のprimitive整数多項式にして計算する。
Sylvester determinantのN!項と各d係数の項数を使う保守的な絶対係数上限であり、
実際の終結式がこの項数・係数高さを持つとの主張ではない。FULLU全域rank3、
空性、full cubic/quartic rankはこのpair選択から従わない。
