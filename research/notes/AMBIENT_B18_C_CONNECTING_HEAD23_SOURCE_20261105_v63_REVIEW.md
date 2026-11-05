# 元 head23 の coherent connecting SOURCE：独立検証

独立な数学 SOURCE 検証は PASS した。報告 SHA256 は
`25cb1e164eb54b85f6835a51781df9cfa0f5adae5279193e22c1a48e34d2b281`。
q7 が読む分子の各係数を sigma0、sigma1_r、sigma1_z、sigma2_rr の畳み込みとして
再構成した。保存済みの二高座標、二 canonical socle 座標、および五 head の低い二行、
計十四個の未評価多項式と一致した。新しい root や出力多項式の値は一つも評価していない。

元の最小 head1 は free56 行23、physical67 の tauA tauB 倍であり、sigma0=e23 である。
sigma1_r=−n0 w の A 曲線行7は −nuB_rr*bB67_rr。
B 曲線は z=−c r² なので行15は −c5(nuA_rr−Σ c_i nuA_zi) となる。
H sector は別の (z7,r²) による表示であり移動させず、行23は delta67_r³ である。
残る H Hadamard 行は delta67_zi*r である。これを全56列へ戻しても、
第一七 F-unit 列の sigma1_r は厳密に零となる。

sigma1_z の全24 source 成分、sigma2 の全24 source 成分を形式変数として保持した。
前者では q 高座標の左 adjoint が Bz0,z4,z5 以外を消し、後者では既受理の L の三行以外を
消すことを独立な畳み込みで検査した。消去前に未知 delta2_zz や source 成分を零とはしていない。
sigma2=−n0(m2 sigma0+m1 sigma1) の両項を保持し、n の一次微分を落としていない。

native V4 sigma0 の A 部分は −nuB2*bB67,2、B 部分は零である。
q の H 二次・zH・B 混合 carry、全 V3 sigma1、原 z=0 の m2 Hadamard carry が含まれる。
従って **この head** に primitive source4 は必要ない。他の head へは拡張しない。
残る新しい label は、既存 joint105 三次級数から合成する bB63 の r ell²、
原 physical root38648 の一つだけである。その値と全出力は未計算である。

free 行19,28,29,72の省略項は、既存 C2 の literal delta30 零と、独立受理済み ALPHA4 の
lambda12 零を追加で pin した。独立 checker の初回は前者を形式変数のまま残し、保存式と
一致しなかった。そのコード・事前 receipt を ATTEMPT1_HISTORY に保存し、実際の零の source
正規化を加えた後に一致を確認した。原資料やその式の訂正ではない。
継続時に出力フォルダの事前 receipt が残っていたため拒否された履歴も保存した。

型の条件付き帰結も確認した。M7=J2²⊕k³ が G4=J2⊕k² に全射する。
単純 source h0,h1 の単純 target k² への像が独立なら、Soc M の像は rG と合わせ三次元以上。
三次元 kernel が半単純ならその像は二次元になるので矛盾し、kernel は J2⊕k となる。
head1 の低二行は零なので、その証拠 minor は low_ellH(h0)*k2(head1) である。
この minor は未評価であり、実際の型をまだ選択していない。

実行は0.070秒、所有上界79,132、係数演算前94 bitで、60秒・2 GiB・150,000 slots・
2048 bit・RAW32 MiB以内だった。値 worker の cache/source/code/custody gate は別の検証である。
本数学 SOURCE PASS だけで実行許可や cache の忠実性を主張しない。
元 H0、K torsion、nil 三乗の結論もない。

再現コードは `code/ambient_B18_C_connecting_head23_source_review_20261105_v63.py`。
凍結出力先への再実行は拒否する。別作業コピーで再検証すること。
