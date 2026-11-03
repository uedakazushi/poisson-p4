# v62 非平方 core：三回目の探索引継ぎ

2026-11-03。旧二つの bundle を保持し、追加作業を `nonsquare_core_wip_v3.tar.gz` と manifest に保存する。これは研究の停止指示ではない。証明済みの重い三次欠損曲線の ambient 除去は既存 v62 theorem packet、J=0 の product 欠損有限性は `S61_BINARY_NONSQUARE_PRODUCT_J_ZERO_20261103_v62` を使う。この bundle の partial pivot、raw modular output、時間切れを全 core の定理と読み替えない。

## 証明済みの追加段階

J=5ad−6e−4=0、adFBC≠0 上では、二変数への厳密な座標表示と四つの QQ 12minor の unit gcd により product-rank41 欠損は零次元以下と証明した。独立 review も保存した。これは J≠0 の開や quartic/ambient の全被覆を含まない。

新しい [斉次 modular certificate 補題](../../notes/S61_BINARY_HOMOGENEOUS_MODULAR_CERTIFICATE_20261103_v62.md) は一般の数学的な検証原理である。各次数の有限 Macaulay map の rank 半連続性と、monic p-integral Gröbner basis の同じ initial module から逆包含を証明する。この補題を実際の大きい matrix に適用する計算は未完了。

## 三 core の現状

J≠0 の mod-ell 24×56 matrix は既知単元で14×46へ縮約される。`thin-r3/bounded-lift.cpp` は一行の column を次数≤Nの単項式で倍し、QQ の有限 Macaulay 係数行列を exact RREF して、既知単元の積を作る polynomial witness を探す。見つかった witness は元の QQ[a,e,d] の恒等式として再検査してから出力する。二項の小 fixture `a*e, a*(1-e*d)` で witness `(d,1)` を確認した。

元14行を次数≤2で走査すると各行に既知単元 aJ 又は aJ² の witness があった。これは各行が局所化後 unit ideal を生成することを示すが、14行全体の独立性を示さない。採用した最初の witness は120項・次数2で、従来 scalar GB lift の199項・次数7より小さい。full-vector native identity と既知単元 pivot により13×46へ進んだ。さらに次の最小行の次数≤3の183項 witness aJ²を exact に検証し、12×45（402839係数項、約256MB）となった。全24 rank の証明は未完了。全行分の独立 witness を先に append しても追加 unit pivot は一つに留まった。

旧 B-free 270×187 matrix には8項の a⁴ witness があるが、貪欲 pivot を進めると約1690万項まで増える。polynomial lift を単に繰り返すことを解決済みとはしない。ZW と ell 商の genuine cubic は元135式から native に再構成した。小さい充分な quartic block の全 core 被覆はまだ得ていない。

## 元135式の全三次 matrix に戻る経路

product 部分だけでなく、現在の135式から再構成した63×91全三次 matrixを使う。凍結 `S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62` と入力SHAを照合する。従来21単元pivotの後の42×70へ exact な定数行・列結合を加えると、adFBCJ開で8個の追加単元pivot（12.3秒）を得た。短い行には次数≤1・50項の J² witness があり、full-vector identity を確認してさらに一個を除去した。累積30pivotで33×62、一般階数の残り目標は27である。一つの厳密な零行を含む。これは全列の被覆ではない。

別に `fresh-cubic/low-syzygies.cpp` は、core係数次数を限定した左syzygyの有限Macaulay係数行列を作る。Fpは独立なequationを選ぶためだけに使い、その部分行列の exact QQ RREF で得た全kernel basisを**元の全列**に掛けて多項式として零か再検査する。この再検査が成功した時だけ出力する。元全matrixのrowを既知単元だけで正規化した後、次数≤2で未知数630・独立equation629、一つのsyzygyを得た。元のrow48の係数は1で、他はrow0,1,2,4,5,6のみである。従ってrow48はこの開の全点で削除できる。

row48を削除した62×91では、次数≤3で未知数1240・独立equation1239となり、別の一つのsyzygyを全係数で検査した。これは二次式の積の間の関係で、元row14の係数12a/5が単元である。従ってrow14も削除し、61×91にできる。階数は57が目標のままで、依存関係はさらに四つ残る。元の63行全体を独立な57行へ全開で縮約した証明ではない。続く次数≤4（2135未知数・80370equation）の native 計算は、120秒・3,400,000KiB virtual cap の範囲内で exit134となり、候補を返さなかった。stderrは有限行列のサイズだけで、abort原因を断定する診断はない。この段階から結論を出さない。

原matrix `full` はheaderの後すぐsparse polynomial bodyが始まる。Schur `schur` には二行のmetadataがある。形式を混同しない。native compile は `c++ -O2 -std=c++17 source.cpp $(pkg-config --cflags --libs flint)`。kernel CLI は `low-syzygies plain-full-matrix degree output [known-units-file]`。`syzygy-fixture` の2×1 `(a,e)` では厳密な `(−e,a)` を確認した。

## 新実行系と未完了計算

既存 Cloud の Singular 共有 library と data を再利用し、同版の Debian console executable/module だけを ignored work prefix に加えた。package URL・SHA、実行 wrapper、既存 library source hash は [実行系記録](SINGULAR_RUNTIME_20261103_v62.md) と manifest にある。既存 setup と network policy を変更していない。

- raw 非斉次 `modStd(M,1)` は J=0 と J≠0 の両方で unit-basis 候補を返したが、公開 library の exactness flag と研究用の逆包含の検証を区別し、全 core 定理に採用していない。
- exact QQ `M*L=E` を要求する modular right-inverse と直接 lift は制限時間内に未完了。
- corrected 斉次 V2 は各600秒まで継続した。J=0 は6,000,000 KiB virtual capで exit14/no more memory、J≠0 は exit124。候補 GB はいずれも得ていない。終了後の全 worker 消滅を確認した。`--cpus=2` でも内部 library は多数の forked task を作る場合があるので実際の cgroup memory を監視する。
- QQ(a)[e,d,z] 上の元14×46の直接stdも120秒で未完了。
- J=0 の四 minor から既知単元を各々除いた二変数多項式は次数69,94,76,69である。M2 QQ GB と直接 Singular slimgb を各120秒まで試したが未完了。四 scalar polynomial の corrected homogeneous modStd も3GB virtual capで候補前に exit14となった。第一多項式は QQ 上既約と native factorization で確認したが、これは四 polynomial の localized ideal が unit である証明ではない。

再生成可能な大きい residual は bundle に重複収録せず、manifest にサイズと SHA を保存する。元 file、recipe、入力 hash を保持している。必要な source/command/log のみを収録し、成功と時間切れの execution receipt を別記する。

## 次に必要な有限証明書

全 marked 非平方を一般 β 上で除くには、実際の bad core の ad≠0 部分を零次元以下にする必要がある。重い三次曲線は既に処理したが、他の genuine cubic/quartic 欠損曲線が残る。product-rank41 の全 QQ 被覆も独立した未完了段階である。

J は c=L=1 正規化を戻すと `Jtilde=5ad−20bc−2ce`、対角 torus character は `(-1,+1)=character(cL)` となる。対応する projective toric GIT quotient では ample section を与える。しかしこの route には実際の bad support の proper な閉包と、Jtilde=0 の全境界点での ambient 除去が必要である。c=0 では Jtilde=5ad、L=0 では Jtilde=5ad−30bc。既存の一般 β 上の c=0/L=0 除去だけを、これら全特殊 core の ambient 冪零性と読み替えてはいけない。
