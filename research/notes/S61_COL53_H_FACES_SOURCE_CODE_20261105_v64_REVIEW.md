# col53 完全 H 二面の SOURCE/code gate：独立検証

独立検証は PASS した。報告 SHA256 は
`7c6e258796e3a5cd5f4b6a5d612c56c2def02b4ee03fd587fbefa11a7ff13783`。
固定 source5e66838…、VALUE gatefe261a…、未実行 worker d7e9a87…の三 hash を pin した。
signed な面係数はレビューで一つも評価していない。

既に独立確認した元 source の1,203 SSA nodeを、今回の H extrema graph と全て照合した。
古い十頂点配列を H の最小・最大に置き換えた部分と、追加の整数 metadata 以外に、
演算・親・literal source・QQ/monomial 変換の変更はない。
raw53 の各再定義を別 SSA とする元 frame を保っている。
キャッシュ cell の九 byte header、353894400000 と shift(1,0,0,0,1) も確認した。

有限 Laurent 支持で積の最小 H は子の最小 H の和となり、その面**全体**は
二つの完全な最小面の積である。最大面も同様である。
和の指定 H 面に寄与するのは H が一致する子だけである。
QQ source factor が厳密に零なら枝を除けるが、既に零と判明した十個の lex 頂点係数を
枝刈りの根拠にはできない。相殺後も source の H extremaを変えず、指定された面を計算する。

この規則で両 SOURCE modeの全五座標支持を別実装で再生成した。
各 node の支持項数・hash、および最終 snapshot の全座標と一致した。
根の H=−3 面は766項、H=21面は214項である。
四座標へ repackせず、符号付き五座標 tupleを最初から最後まで保つ。

両 VALUE modeの依存集合、H-tied引数、全 last-useを独立に再構成した。
source q=0以外の追加枝刈りはない。全876個の整数 rescale 恒等式と、
実 VALUE経路の388個の convolution係数 majorant を確認した。
copy・積・和を固定 source q倍の整数 mapとして計算するため、
各項と部分和の大きさは source M と整数 weightで事前に抑えられる。

decoded D/P/F は28,362/6,681/4,904で、numeric label、index、必要集合、
allowed tuple、ROOT document、failure予備を含む META は98,171。
全五座標・係数・索引に一項7 slotsを課し、最大 live+workは2,912項・26 maps。
leaf decode16,152を同時に含む実 code の将来正値上界は126,262で、hard150,000内である。
mapをcacheへ移す際は `cache[k]=work; work={}` とし、複製しない。

通常 committed640 bit、temporary1921 bit、hard2048は変更していない。
整数演算の source一時上界と、codec の保守的な naive分子上界1372 bitを確認した。
codecは全五座標・primitive gcd・有理数contentを streamで往復し、空 mapの九 byte零表現も扱う。
出力の完全な cell を一時ファイルから renameし、途中の mapを完全な面とは主張しない。

55秒 signal または失敗時は live cacheと部分workを stream保存する。
既存 map はコピーせず、failure予備8192を一時的に metadata chargeから外して、
alias key listとcodec領域の最大6906を確保する。
indexは配列でなくJSONL stream、生成物2 MiB、全fresh RAW17,661,238 bytes／32 MiB内である。
math50秒・signal55秒・hard60秒・2 GBの範囲と、自動 retryしないことを確認した。

独立 reviewer自身の SOURCE支持計算は0.593秒、最大所有見積94,587であった。
この PASS が親の既存条件付き **同じ一回の二面 job** の条件を満たした旨を担当へ通知した。
実際の面の値は、実行後の別 VALUE peerを要する。

両面が零でも成分全体の零は証明されない。非零の面は固定成分の非零を証明するが、
全167×14・87追加RHS・66制約・oneP290・F9/F10・全補集合についての
rank、inverse、分類の完了はこの gateから主張しない。
parameter pointやprimeを追加していない。

再現コードは `code/s61_col53_H_faces_source_gate_review_20261105_v64.py`。
凍結出力先を変更せず、別作業コピーで再検証すること。
