# v63 CK：元接続写像・T2指定二次関係・一般五点の非零entry

2026年11月6日（日本時間）。直前の通常push済みmainは`1004a1f7212dc0ae53a063a1effe4a161fa9bbda`。Cloud開始SHAは`4c3f623fb698f9adc6077015b25e2325b11072af`。

以下の八項目は、それぞれ明記した加群・有限商・開部分に限る厳密な証明または厳密計算である。関数体上の非零係数と大域的な尽くし、有限体gcdと標数零factor、比較Cと元H0を区別する。新規性は独立調査未完了。

## 元接続写像の核J2+kと七次元商での分裂

元q7全列と全mixed carryを残し、既存105変数の三次以下cacheだけからhead23の全十四係数を求めた。lambda=2039/360、正規化highは(-55495/968256,1331/9216)、低行はh0=(4/5,2/25)、h1,h2,h3=0、h4=(-9396/42025,-4248/210125)、chi4=-36/1681である。別のbinomial座標法による全十四値の独立再計算が一致した。単純商の二次minor1331/115200は非零なので、M7からgenuine G4への三次元核はJ2(r)+kとなる。R0安定なS=span(h4,rh4,h0,h1)の像の行列式11979/2260608800は非零で、未評価の八high parameterに依存せずM7=S directsum (J2+k)を得る。

- 分裂はM7=Told/(z7,r²)Toldだけ。Toldへのsectionや元H0の同定を含まない。
- 残り八high parameter、Kbarのextension、元H0・nil三乗・指数は未決。

証拠：

- [`reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHE_SOURCE_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHE_SOURCE_20261105_v63_REVIEW/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63_REVIEW/report.json)

## 元Told接続写像のR0非分裂

Toldの実際の次数二初期形でr²はh4の定数係数だけを独立なq2=r²h4へ送る。genuine G4の巡回生成元のliftはh4定数係数1/chi4=-1681/36を必ず持つ。高次補正のr²像はm³に入るので、このq2係数を消せない。従って元のTold長十一からG4への全射にはR0加群sectionが存在しない。独立filtered-module証明と全係数点検が一致する。

- 七次元商M7での分裂と両立する。Kbarの全作用表や元H0の長さは結論しない。

証拠：

- [`reports/AMBIENT_B18_C_CONNECTING_TOLD_NONSPLIT_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_TOLD_NONSPLIT_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_TOLD_NONSPLIT_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_TOLD_NONSPLIT_20261105_v63_REVIEW/report.json)

## 接続写像によるr三乗の因子化とKbarの消滅

元Toldの五head全てのr³像を再構築すると、h4だけが非零deep classを与える。qC=z5*r²h4=0等の元関係からr² ITold=0である。従ってr³:Told→k*gdeepは(gdeep/chi4)*chi*deltaを経由し、接続核Kbarではr³Kbar=0となる。独立実装が五像と全元関係の使用を確認した。

- 元Kへ戻せる結論はr³K subset ell*Kだけ。r³K=0、H0=0、nil三乗零は結論しない。

証拠：

- [`reports/AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63_REVIEW/report.json)

## native C28のr二乗ell三次係数全34個

元raw33からP7への全1089積と同じcanonical枝を保ち、free-four trace拡大C28で八G8のr²ell全32成分と二monicを計算した。重みの最大五からradical六乗零となり、P inverseは三項で正確である。算術前37bit/temporary112bit・owned14403の証明と独立SOURCE点検を経て、一回計算の全34値を元C132のrational rootsと別のtrace基底積表から独立に再構築した。二monic値は共に零である。

- 同じcacheからのunit関数値だけ。native Xp/Xc/N三次やT2 upper像は別であり未評価。

証拠：

- [`reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_SOURCE_V2_20261105_v63/report.json`](../../reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_SOURCE_V2_20261105_v63/report.json)
- [`reports/AMBIENT_B18_NATIVE_C28_ROOT_SOURCE_PEER_20261105_v63/report.json`](../../reports/AMBIENT_B18_NATIVE_C28_ROOT_SOURCE_PEER_20261105_v63/report.json)
- [`reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_VALUES_20261105_v63/report.json`](../../reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_VALUES_20261105_v63/report.json)
- [`reports/AMBIENT_B18_NATIVE_C28_DIRECT_REVIEW_20261105_v63/report.json`](../../reports/AMBIENT_B18_NATIVE_C28_DIRECT_REVIEW_20261105_v63/report.json)

## T2全617 upper式のOLD二次形式係数

全757 typed originと全617 kernel recipeから、完全なGamma(Xc3-Xp2*a1)+D(Xc2-Xp1*a1)のOLD二次線形感度を構成した。1926候補slotのうち1909が形式的に非零、106case-channelに42functionIDが現れる。別実装は元757 differentialをDで先に収縮し、その後617 source和とliteral outer補正を行って同じ符号付き1909係数を得た。opaque第二記号に依存するのは十二caseであり、残り605はこれら十二記号から独立である。

- 605 caseが零または計算済みとはしない。gA/gBのmixed次数三certificate、physical秒次、known-FIRST cross、native carryを保持する。
- 十二opaque関数は線形第二感度の話。D multiproduct祖先やknown-FIRST値から消えたとはしない。T2 HFと全長は未決。

証拠：

- [`reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63_REVIEW_V2/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63_REVIEW_V2/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SCOPE_CLARIFICATION_20261105_v63/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SCOPE_CLARIFICATION_20261105_v63/report.json)

## T2の指定upper508と288が厳密に零

case508のOLD595単独sourceはa1とnative c0が零で、coherentなB21 fiber-product商により元全九firstと必要r²ellを保つ。case288はOLD352単独、a1全八entryがz1のみで、純plane Gamma Xp2*a1は零だがD_z1r Xp_rとD_z1ell Xp_ellの既知FIRSTcarryを全て残す。それぞれGamma三次、D二次、case288全八FIRSTcarryの非負係数envelopeが両channelで個別に零となることを独立dense56実装で確認した。従って符号付き値workerを走らせずに、この二つの実upper classが零と証明される。

- 指定二caseだけ。残り615 upper、native549、全605case、T2 HFや長さを結論しない。

証拠：

- [`reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63_REVIEW/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER508_INPUT_GATE_20261105_v63/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER508_INPUT_GATE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER508_INPUT_GATE_20261105_v63_REVIEW_V2/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER508_INPUT_GATE_20261105_v63_REVIEW_V2/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_20261105_v63/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER288_INPUT_GATE_20261105_v63/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER288_INPUT_GATE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_INPUT_20261105_v63_REVIEW_V2/report.json`](../../reports/AMBIENT_B18_ELL2_T2_UPPER288_SOURCE_INPUT_20261105_v63_REVIEW_V2/report.json)

## 一般五点col53固定entryの非零証明

元全leaf decoderと1203 SSAから完全なH=-3面とH=21面を全五座標を保持して求めた。sourceのq零とH極値が一致しない和だけを除外し、既知vertex値が零という理由では除外しない。H=-3面の実非零係数は411個、H=21面は全て零である。指数(-3,0,6,-1,10)の係数241815/134217728が非零証人となる。別の逆順整数replayが全411係数と零面を再構築した。従って元full F11は非零で、既存の非零prefactorを用いたp=P7*F11/(U²*F10)も関数体上で非零である。

- 411項は面だけでありfull F11ではない。十一番目typed/元83番目中七十番目pivotの存在と将来shape166x13はSOURCE帰結だけ。更新済み行列・RHS・inverseは未計算。
- 全F9/F10/F11補集合、87追加RHSと66条件、P290一回を保持。一般五点の尽くしや全分類上界は未決。

証拠：

- [`reports/S61_COL53_H_FACES_SOURCE_20261105_v64/report.json`](../../reports/S61_COL53_H_FACES_SOURCE_20261105_v64/report.json)
- [`reports/S61_COL53_H_FACES_SOURCE_CODE_20261105_v64_REVIEW/report.json`](../../reports/S61_COL53_H_FACES_SOURCE_CODE_20261105_v64_REVIEW/report.json)
- [`reports/S61_COL53_H_FACES_VALUE_GATE_20261105_v64/report.json`](../../reports/S61_COL53_H_FACES_VALUE_GATE_20261105_v64/report.json)
- [`reports/S61_COL53_H_FACES_VALUES_20261105_v64/report.json`](../../reports/S61_COL53_H_FACES_VALUES_20261105_v64/report.json)
- [`reports/S61_COL53_H_FACES_VALUES_20261105_v64_REVIEW/report.json`](../../reports/S61_COL53_H_FACES_VALUES_20261105_v64_REVIEW/report.json)
- [`reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64/report.json`](../../reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64/report.json)
- [`reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64_REVIEW/report.json`](../../reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64_REVIEW/report.json)

## 固定有限体gcdと標数零共通因子のa次数上界

固定p=32003、a/e/dのKronecker重み646/19/1で二入力のgcd次数は3899、monic支持は94項となる。独立検証は元二入力をこの候補で全係数割り、余り零と二商のEuclid互いに素を確認したので、候補は厳密にこの有限体gcdである。元整数入力のleading係数がpで非零のため、Gauss補題による任意の標数零共通因子のKronecker leadingもpで生き残り、次数は3899以下である。因子のe,d座標次数箱内では重み646が支配するので、共通因子のa次数は六以下、a次数六ならleading e,dに19e+d<=23が必要である。

- 非定数な標数零共通因子の存在・互いに素・radical/Fitting・補集合消尽は未証明。CRT四十素数のSOURCEは準備であり値計算未実行。

証拠：

- [`reports/S61_TWO_TAIL_FIXED_KRONECKER_SOURCE_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_FIXED_KRONECKER_SOURCE_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63/report.json)

## 保存した未完部分と次の具体的操作

native二十三pivotの全九十二Gamma Xp二次値は一回計算済み（primary `8eea9b07…`）だが、異なるuniversal trace法による独立値確認はこのcheckpoint時点で未完了。全L二次S零次とmonic carryを保持する。case549はOLD298/4・native649の和で、native649のS零次零を正確なinputmaskで証明してから全functionalを評価する。全617の上位像やT2長17〜20は依然未決。

HEAD29では元delta63=g20144−g40922の未知105変数に対するorigin gradientが全て零と、reverse ADと独立forward first-jetで確認できた。このため選択四次deltaは既知三次以下Phiだけに依存する。HEAD30のmixed/zzにも同じ原理を適用できるか、有限六十三次元代数の安全出力・lower extensionをSOURCEで証明中。未計算Phi四次などを真の零と仮定しない。未評価のhead高位八parameter、Kbar extensionと元H0は残る。

col53の十一番目typed pivotはfull F11 circuitの非零から存在するが、新166x13係数行列・RHSやinverseは未計算。正確なfactored callbackと全所有量/高さ/RAWを先に設計する。旧十vertexの値は独立再計算済み。旧producerの84916という予算表示にはtuple label40100が欠けていたため、訂正receipt `S61_COL53_VERTEX_PRODUCER_BUDGET_CLARIFICATION_20261105_v64`で因果的なdecoded VALUE上界125016を別記し、凍結済みprimaryを改変しなかった。H面の新workerでは全所有124092が算術前126262以内、独立peerは128188が独立上界130358以内である。

二tailの固定四十素数CRTはSOURCEのみを保存し、候補をprimitive化した後の条件付き8192bit境界と二本の整数割り算・全係数掛け戻しを準備する。有限体gcdが正次数でも標数零因子はまだ得ていない。

元nil指数3〜6、元H0の長さ・ell消滅・nil三乗、T2全長17〜20、全FULLU57/239、全補集合、一般五点の尽くし、全分類の上界は未解決。下界58（rank4は52）、NEW式28次と五つの凍結原稿は維持する。続行可能な研究を保存後も継続する。
