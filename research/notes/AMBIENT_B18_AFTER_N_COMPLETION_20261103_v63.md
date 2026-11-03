# B18 frontier引継ぎ：N全完成環は解決、CA/CBはnu欠損とcubic classまで

2026-11-03。凍結した旧14ページ稿と全台帳は本agentでは変更していない。
親agentが証明status・索引・Git保存を管理する。

## 解決したN

actual N odd2jetの独立reviewから、既知J=(I,nu)の全9 minimal classをactual Iが張る。
NakayamaでI=Jが確定した。

    R=C[[g19,nu,s6,w,v]], I=(nu,s6,wv,v²),
    A=C[[g19,w,v,s6,t²,s6t,t³,wvt,v²t]].

cone dim28/edim37/depth20/e2、projective dim27/edim36/depth19。
全81関係（36 Taylor syzygies＋45 products）で全kernelを与え、別M2 kernel計算もPASS。
Aはdomainであり、embedded primeはconductor quotient R/Iにのみ現れる。
t²,t³∈A、t∉Aなのでnormalizationはelementary subintegral/universal homeomorphism。
どの次数でもoffdiagonal geometric/reduced-arc pairはない。
全N-germ、source critical scheme、commonR29へのpartial pinch projectionを区別する。

* `AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63` と独立REVIEW
* `AMBIENT_B18_N_FULL_COMPLETION_20261103_v63` と独立REVIEW

以前の八生成N ansatz二つは、actual minimal生成数9で除外された。
full paired-flag三次obstructionを凍結する必要はなくなった。
ignored workのfull9flag計算はNのsubintegrality定理に使っていない。

## CA/CBで確定した違い

actual全odd一次rankは7、discriminant tangentへの制限rankは6。
rank6はadjusted eta一定であり、raw row61固定のraw Y28列はrank7である。
source Jacobian全150を使ったactual二次項では、zero-R-firstjet generatorsのspanは

    span(WU,WV)

のみで、UV classが欠ける。
一方chosen nu-liftはCAで−UV/75、CBで−16UV/9を持つ。
従ってnu∈Iならi_nu−nuのzero-firstjet UVが矛盾するため、**nu∉I_CA/I_CB**が独立受理済み。
Nと同じI=JへのNakayamaは成立しない。

* `AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63` と独立REVIEW

## 新しいcubic witness（author PASS、独立reviewへ引継ぎ）

全七monic critical消去のU-axisを使うと、全150 vertical derivativeは次数2まで零、
branch61も三次まで零となる。
七つのactual linear基底と二つのactual quadratic基底をsubtractした残りはraw row63だけで、

    CA: (12/125) U³,
    CB: (−4/243) U³.

actual IのJ/mJへの像は7+2+1=10次元、J/mJは11次元である。
mu(I)≥10だが、I∩mJをmIと同一視して上界10を結論しない。

* `code/ambient_B18_CA_CB_odd_cubic_class_20261103_v63.py`
* `reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63/`

## 撤回した方法と残る型

R19\NのCA旗とCB旗の相違を、個々のnormalizationの二preimageと混同してはならない。
既証では個々のCA/CBはR19\N上で滑らかである。
このR19曲線を個々のconductorの追加枝として使うmultiplicity下界は撤回され、台帳にも採用していない。

leading dataの候補は

    I=(s6,nu+cUV,WU,WV,U³)

だが、full Iの生成上界・nuのnilpotence・exact depth・全次数formal同値は未解決である。
追加の高次generatorがV-axisの厚みを減らす可能性を排除できていない。
actual D=CA∩CBの全完成型・primary構造も依然未解決である。

## 次の有限な計算候補

critical V-axis U=W=0の正確なsource係数は、CAでは

    (s0,...,s9)=(V/10,0,0,V/10,−V²/200,0,0,−V/20,V²/400,−V³/24000),

CBでは

    (V,0,0,V/2,−V²/4,0,0,−V/2,V²/4,−V³/24).

V=1でのGL orbit rankはCA20、CB19（いずれもscalarを含む）だった。
これは有限QQ点の測定であり、CBの一般V orbit次元を外挿しない。
sourceの正確なphysical点は`work/cloud-manuscript-20261103/CA_CB_Vaxis_point1.json`に保存した。

generic V-axisでのnormalizer kernelに沿うactual odd初出次数p³,p⁵,p⁷を、
全source方向を保ってcommon target座標をimplicit固定して計算する方法が候補である。
候補ならV-unit枝はcusp(2,7)となりconductor長3、W-unit枝はcusp(2,3)で長1になる。
この長さ4のlowerと、三次initial上界4のCM性が一致すれば、有限module比較へ進める可能性がある。
現在はこのgeneric cusp-order・長さ下界を証明していない。
