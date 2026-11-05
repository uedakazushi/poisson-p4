# HEAD29/30 joint V2：混合精度 SOURCE・コード・全資源監査

独立監査は PASS した。報告 SHA256 は
`985e63285949c4503136288971c13494f898ae15e6e97325b0d17d5b87e6ec5c`。
SOURCE fabde8dd…、compiled 2d8b4541…、profile aff1e8a6…、
worker faf2c7f4…、custody 7e64f0bd… を固定し、signed出力は一つも評価していない。

元のscalar記録を別readerでたどり、全8081 nodeと全child、root、定数n/dを照合した。
delta二rootの祖先4387個だけ63精度、残る3694個はplane10である。
63 nodeが10 childへ依存しないこと、全retirement順序、最大241 live／8662 QQ slotを再構成した。
元105 lower入力は受理済みB63の既知31／未知socle32の拡張と文字通り一致する。
未知socleのzero skeletonは真のimplicit係数が零という意味ではない。

内部で63精度を要する共有bB63 root38648も、bB出力としてcaptureする時点で
**63→plane10射影**することを、sourceとworkerの両方で確認した。
delta差全63だけがsafe出力であり、bB八rootはplane10だけを採用する。
旧V1のprecision capture失敗は保存し、V2と混同しない。

元producerの係数逆元recurrenceと違い、独立監査はnilpotent正幾何級数で
逆元majorantを作った。全8081 nodeの正majorant hashが一致した。
両逆元積diagnostic、原点・raw n/d parse、lower overlap scalarを含むroot上限51／pre103を確認した。
全八未評価postoperatorのmajorantと、追加canonical checkのpost118／pre237も確認した。
head30のBは118 bit、scaled Aは113 bit以内なので加算前237 bit以内である。
約分後の結果は係数ごとに証明したcanonical operatorそのもので、その別envelopeは116 bitである。
初回の監査側はこの相関を捨てた独立絶対値和119 bitを課して止まった。
そのcodeとdiagnosticを保存し、正しいoperator相関を明示して完了した。
worker・SOURCE・private capを変更していない。

numeric labelsも含むruntime1は67465で、元記録・精度・retirement40184、
105入力6720、live QQとkeys9144、rootsと積table1216、lowercheck201、予備10000からなる。
runtime Counterは無く、input配列はin-place変換する。
phase1を保存してから元record/input/product bufferを解放し、postoperatorを読む。
postの全2673積項・RAW生成物・guard時記録を含めfresh RAW上界1064308 bytesである。
standard-libraryのみのworkerをASTで読み、無関係な有理数guard fixtureも独立確認した。

独立SOURCE phase0上界110195、majorant算術上界97716、実行1.964秒だった。
実行する場合は親の既存限定ONE GOの下で、正確なfaf2 workerと当該peerを指定し、
`prlimit --as=2000000000 --cpu=60:60` を用いる。
384/1153・global2048、hard100k/global150k、RAW2 MB/global32 MiBを維持する。
phase1が完成した後にpostが止まっても、それを保存し自動再実行しない。

この監査は準備条件であり、actual root/head値の検証は実行後の別peerで行う。
全Kbar作用、extension、元H0、nilについて新しい結論を与えない。
再現コードは `code/ambient_B18_C_connecting_heads2930_cache_source_review_20261105_v63.py`。
