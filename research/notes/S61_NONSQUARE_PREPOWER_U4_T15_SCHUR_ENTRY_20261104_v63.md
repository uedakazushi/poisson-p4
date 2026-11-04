# direct Schur 一成分 u4・target15 の実計算

同一の型付き768演算と精密化したsource資源証明に従い、compact行4/current15行5、元32標的15の cleared numerator 一つだけを評価する。Mark の HYBRID residual6×13 位置(0,9)との関係は

    residual_(0,9)=rho_u/(delta2*rhoLCM*D4)*Saux

である。source-binding.json は全原子、演算、同一pair-only H行とP4恒等式からのこの等式を固定する。h22 の逆元も下位H3の逆元もこの一成分にはない。全六行の chart FULLU∩D(Delta5 R434) とその補集合は保持する。

全入力原子は固定されたQQ係数とaのLaurent冪だけである。今回の全768ノードのB分母冪0、未計算原子のliteral分母形状も検査する。他の連成列にある99個のB分母は今回のclosureに存在しない。係数計算を省略して可逆とみなすことはしない。

入力保持を含む予測は143283+5895=149178項、上限150000。native は実際の全生存項、literal reserve5895、一般reserve15000、output-size scratchを毎演算で数える。aの指数をin-placeで動かすと同時に外のa冪を逆に動かし、各項の実指数保存を確認する。二つの入力の丸ごとshiftコピーを作らない。

完全な演算DAGと全原子が再現可能なsource convolution恒等式である。各768ノードの項数とrepresentation shiftを逐次記録し、三つのH源cofactorと最終Saux・cleared numeratorの全QQ係数を保存する。巨大な全中間多項式の複製を保存する指示ではなく、独立検査は固定768演算を直接再構築できる。

この成分は零又は非単元にもなり得る。非零でも全六行の階数、全九核、global full57は従わない。別の成分や小行列式、点、素数、一般因数分解へ自動的には進まない。結果は同名reportに固定し、ガード失敗では部分資料を保持して停止する。

実行上限60秒・2GB・32MiB・150000項。元のBOX失敗資料と精密化preflightは変更しない。
