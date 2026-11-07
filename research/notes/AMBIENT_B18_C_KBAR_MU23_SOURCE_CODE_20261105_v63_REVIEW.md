# 元 K の二生成・三生成判定：SOURCE と V2 実行前検証

既存キャッシュだけを用いる37本の二次係数と十四本の判定式について、SOURCE・枠・コード・保存方法を独立に照合した。requested rr 係数の抽出、37本の値との内積、十四作用の値、非零 minor の値は実行していない。

既知の rank(r|Kbar)=4 により Kbar/rKbar は U,S,q2 の三次元である。十四本 z_iU,z_iS の一次 NF は実際にすべて零で、商では q2 方向だけが残り得る。z_iq2 は旧 Told の m³=QQg に入り、g=rq2 なので商で零となる。従って十四係数のどれかが非零なら最小 R0 生成数は二、全部零なら三。ell を法とした有限生成 module の Nakayama により、これは同じ元 K の R0 上の最小生成数を判定する。D 上の七生成、元 H0、nil 三乗の判定とは異なる。

元の選択37本と正規化した37本の間の全1480成分を multiplyback で照合した。十四の一次式、および rU,rS は、producer の単項式除算を使わず、37×37の線形方程式を下側 pivot 順で消去して求めた。定数 lift、825列中の168+657という配置、元768列の実際の constant combination と同じである。37本それぞれの full56 X0 は全2072成分で零。

完全な q2 dual は literal index179、すなわち旧 head4 の r² 係数である。全40個の初期 relation から生じる二次関係をこの dual が消すこと、既存の exact HF2=2 と標準二元 q1,q2 による完全性を用いる。単なる斉次 relation を真の M2 に置き換えていない。

直接 section の元55行は全次数で零、Pi0 の最終行は元55行そのもの。全640 raw first、選択27 pivot first の55行 r 係数、84本の unit-cycle の head4 r 係数もそれぞれ零である。したがって first/second Schur carry を保持して展開しても、必要な真の M2 は元 C0 の P9 最終行の二次係数に帰着する。選択37本の全56 X0=0 も直接確認し、未知の二次 projection と非零 X0 の積を落とす近道を使わない。

必要なのは元 P9 の21個の literal scalar root、合計98 source weight である。この21 root と十四本の二段合成・直接合成の重みを元 source から独立に再構成して、compiled317b…と一致させた。提供済みの多項式全体と opaque-second gradient が空であることを比較し、同じ actual C0 の hash と元8変数の monomial 枠を確認した。rr の一成分と重みの pairing はしていない。既存の非線形 CB transport を消したり、未提供の高次係数を零としたりしていない。

SOURCE 検証中、inventory の二行、z5U と z5S の degree2_explicit_source が空配列になっていることを検出した。producer の reducer が、既に標準形の q1 項を入力辞書から pop したためである。真の F2 は −cU*q1、−cS*q1 で、保存済み degree2_known_NF はその非零値を正しく保持する。q2 成分はどちらも零なので、十四の compiled functional と将来の判定値は変わらない。元ファイルを変更せず、この二行だけを別の補正証票に結び付ける。独立 reviewer の最初の SOURCE assert とコードも draft history に保存した。

全 provided coefficient に対する正の numerator／denominator envelope を使い、37本、十四の直接合成、十四の二段合成、lambdaS² を掛ける witness minor の上界を再構成した。将来値は committed136／temporary273 bits に収まる。実際に SOURCE で扱った値の最大は stored67／pre119 bits。

V2 worker94d50…は compiled の hash を固定し、37本の origin・r first の零照合を rr 選択より先に実行する。phase1 を保存してから十四本に進み、二段と直接の二通りを全項照合する。SCI は parser、lookup、積、加算の pending 操作、直前 accumulator、途中の term、computed candidate を保存する。wall timeout を含む例外でも scientific prefix を残し、自動 retry・cap 引上げを行わない。旧 V1 は未実行のまま保存する。

独立な shape-only fixture では未知値を無関係の136-bit文字列に置き換えた。SCI は2108数値、保守的な save 全体8555、failure 全体10202で、いずれも固定 static14076／hard60000より小さい。後二つは独立 reviewer の追加保存・controller 余白も重ねた値であり、producer の因果順序の見積りを変更するものではない。10,000文字の traceback を含む prefix は143,824 bytes、whole510,000／hard800,000 bytes 内である。実キャッシュ値を用いない candidate guard fixture と --help も確認した。

独立 SOURCE 所有上界は68,724数値、全入力と出力余白は14,237,051 bytes。60秒、2GB、global150,000数値、32MiB の範囲内である。将来の値計算は親の ONE GO に従う別段階であり、このレビューは元 K の生成数の実際の二・三いずれかをまだ主張しない。
