# 全605関係の complete 正値 SOURCE 上界の独立査読

primary は `UPPER605_POSITIVE_WHOLE`、SHA-256 `9431e686c987488602817964dcaa3bee47177a2684c776314478a16dfea45a8c`。独立報告 `reports/AMBIENT_B18_ELL2_T2_UPPER605_POSITIVE_WHOLE_20261105_v63_REVIEW/report.json` の SHA-256 は `f204d77f1ccf2fcd1d7ae78449a5a8dc47ab9ac7493229e04931ef2988f6fab4` で、**complete 正値 SOURCE 査読は PASS**。この検算は非負上界だけを計算し、605関係や残る20関係の符号付き応答値は計算していない。

最初に cache/source の同一性を独立に調べた。96境界のうち48関数は literal physical164 の raw33 値を実 B21 に射影したもの、残る48関数は実 FULL9 一次係数と任意の高位拡張である。全9方向の一次を、元の独立 full9 source カタログと照合した。後者の高位零は真の jet の零という意味ではなく、全形式消去証明 `36a453a1…` を適用した**訂正後605関係だけ**で使える拡張である。gA/gB の真の三次を一般に零とする論法は使わない。

ordinary B21 について、raw33 の全1089積と全660 NF が21 monomial への射影を保つことを今回も独立に検査した。FIRST は元 old640/native128 の full9一次 cache と現在の27 pivot に結び付けた。全162個の列の r・z 一次係数、native ℓ 列の `L0 Se+Lℓ S0` の両項と全個別 S0 を保持した。全619個の source 組合せ、full `a1`、563/609を引く係数を元 source と一致させた。受理済み native92 の全23 source と rr/rℓ の両 channel も、独立 ROOT92 査読の実際の値へ結び付けた。

正値再計算では、係数番号を添字とする生産者と異なり、9変数の monomial tuple をキーとした。各和は逆順に経路を集め、未約分分母の LCM を求めた後で絶対値分子を共通分母へ移した。この方法で literal C0 532ノード、全 typed 757ノード、全324 FIRST表、raw619の四部分上界、563/609訂正後605関係、native72 FINE、native92外側、Q の全共通上界を厳密に再構成した。D/U の商だけでなく分母証明そのものが一致した。

boundary/C0/typed/FIRSTの正値係数セル数は152/119/324/1295で、Boolean fingerprint と全324 FIRST表の項数にも一致した。ordinary rr は一回であり、多重積からの既知 FIRST 交差項を落としていない。raw619の部分上界は、563/609を含め形式的拡張の部分値として保持した。訂正基底の upper を零と置いていない。

OLD4 の削除は、受理済みの exact12関数線形性と全形式係数消去が**訂正後605関係**で適用できるために限る。raw609の carry の不在を仮定していない。native75の補助値は形式的、kernel72制限後だけが実 functional という `0384d262…` の範囲を保った。native の coarse354を直接 Q に掛けず、各 FINE 上界へ戻してから whole upper と Q を組み立てた。Q は Ψ・full comparison 補正を含む受理済み source `d6982e39…` の同じ行列である。

605関係を全項込みで調べると、正値上界が非零となるラベルは正確に次の21個だった。

`525,529,530,550,552,554,558,560,562,566,567,569,573,576,577,589,595,599,603,612,617`

それ以外の584関係では complete 上界が `(1,0)` なので、真の二つの upper と Q による六つの class が零であることを SOURCE から結論できる。各605関係についても六つの Q 正値上界を作り、共通上界による分母可除性と分子支配を検査した。上の21個は**実非零という意味ではない**。554は別の受理済み独立 VALUE 査読で実際に零と確認済みなので、残る20個は別の固定 VALUE ゲートと ROOT の許可を待つ。

whole 正値 SOURCE の演算結果最大257ビット、演算前の一時チェック258ビットだった。native 中間を含めた committed354/temporary1063を維持し、hard384/1153を上げていない。演算前 whole owned 上界127,392、観測保持上界96,500、実行0.346秒で、150,000数値、32MiB、60秒、2GB内で完了した。入力 stage は廃棄時期を分け、C0 jet をその表へ縮めてから捨て、typed jet も FIRST/native へ進む前に捨てた。

この査読は584の構造的零と complete 正値上界を閉じるが、残る20応答の値、階数・HF・全長・元H0・nil の新しい結論は与えない。opaque12は605の外に保持する。whole VALUE worker、custody、停止時 codec と寿命、ROOT の新しい実行許可は別の要件である。primary・既存コード・報告のバイト列は変更していない。
