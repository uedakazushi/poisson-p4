# R434 を反転しない direct P5＋H5 の5×12表示

2026-11-04。新 scalar coefficient/minor/rank は計算しない。先行第五pivot packet と
直接P4＋H5、指定一成分の受理済み768 source演算を入力とする。前packetは凍結保持する。

## P5は実 source 行から作る

元15×22をMとし、u=current15 row5、pair rows2,4とする。既存pair-only source演算の
三係数を c_u,c_1,c_3 と書く。その逐語定義は

    c_u=rhoLCM*delta2/rho_u,
    c_1=-rhoLCM*mu_0/rho_1, c_3=-rhoLCM*mu_1/rho_3,
    (mu_0,mu_1)=N[u,pair]*Adj(pair2), rhoLCM=a*F*C²

である。ここでmuとrhoを別規約に置換せず、先行768演算の実nodeをそのまま参照する。
row5を `c_u*M5+c_1*M2+c_3*M4` にする。c_uはD(Delta5)上の既知単元であり、
row2,row4は変更しないため逆操作も明示的である。新row5のheavy五成分は全て零。
残るpure17成分は、先行一成分演算に現れたBhat行そのものである。

次いで従来と同じ五つのH sourceだけrhoで割り、heavy target31をgammaで割る。
P5 source順はcurrent15 `[0,1,3,6,5]`、pure target順は元32
`[22,11,13,25,15]`。P4は従来と同じ四行四列である。

    det(P5)=D4*Bhat[15]-Bhat[piv4]*Adj(P4)*P4[15]=Saux.

これは普遍的 bordered determinant identity であり、D4≠0でのSchurから多項式環へ
戻せる。従ってD4又はR434の逆元を必要としない。Sauxの7022係数は既存primary、
非零性は独立symbolic32003 peerから参照し、新たな5det展開を行わない。

## 指定された単位領域と完全operator

H5 sourcesはcurrent15 `[2,4,7,13,10]`、heavy順は元32 `[27,31,28,30,29]`。
新pivot10は

    K10=[[P5,0],[C55,H5]], det(K10)=Saux*Delta5.

既存のH5 inverseはpair2＋B3のcofactor recipeであり、依存closureを調べても
R434,D4,h22の逆は含まれない。P5inverseは未展開Adj5(P5)/Sauxを使う。
新しい領域は厳密に `FULLU∩D(Delta5*Saux)`。R434=0を除外しない。ただし
R434=0側に新しい非空点を示したとは主張しない。

全row replacement、hybrid unit、二permutation、Schur左右演算と各逆を保存する。
元15×22に作用させるとdiag(I10,T5×12)となる。残る5sourceと12target、全8Q9tailは
先行第五pivotと同じである。元25sourceへの輸送もR434を反転しない既存prefixへ戻す。
新しいscalar unknownやq-dependent inverseはない。

元FULLU sourceのfullrank十分条件は、指定領域でI5(T)=1と同値である。全5rankは
未計算。Saux=0、Delta5=0、p=0の他chartは未処理であり、H7全域をこの一chartで
覆ったとはしない。特にP5の非零detからglobal cubic57は結論しない。

## 次の具体的な計算前見積

先行source支持だけの最小bordered6候補も、output上界121500＋入力58979=180479で
150000を超える。したがって新scalar jobは提案しない。7022項同士を掛けるより、
直接bordered6、又は二つのpair-only行を元のheavy2列へ戻した8×8式を使う方法を
保持する。次のscalar値には、別のlast-use/degree/height証明と明示予算が要る。

R434の一般factor/gcdはまだ試していない。全pでの厳密除法は非零remainderで止まり、
proper p-factor共有を排除していない。このpacketではRの実434係数からdegree、
exponent box、QQ content clearingだけを読むpreflightを保存し、factor/gcdを実行しない。
R0の解決はP5chartの局所化条件を弱めることとは別の問題である。
