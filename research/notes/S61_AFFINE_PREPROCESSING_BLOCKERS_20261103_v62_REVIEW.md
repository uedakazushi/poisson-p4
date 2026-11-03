# 全134二次係数と有界実行の制約：独立監査

2026年11月3日。pencil producer と既存失敗receiptを保存するだけのproducerを、別directoryで再実行した。
失敗した全GB自体は再実行していない。数学的成功とreceipt保存成功を区別した範囲を受理する。

独立コードは元の保存134式から C(s)=C0+ΣsiCi の全有理二次係数を再構成し、保存modp配列と一致させた。
特に **QQ上の** Ci非零行数が各64であることを確認した。
modp非零行数だけからQQ rankの上界を推測せず、これと非零64minorによって正確なrank64を得る。
25個のparameter-free列の有理span恒等式と16次非零QQminorを検査し、正確なrank16を受理した。

選択134minorのmodp逆行列と五つのNi、非零traceを独立に計算した。
M0のdetはp単元なのでQQ逆行列もp整、Jacobi determinant微分公式から各係数方向でdetは非定数である。
従って **この選択minor** の全Niをstrict triangularとして有限nilpotent逆和を使う近道は不可能。
他のconstant-determinant minorや別の前処理を排除した結果ではない。
画像flagの134→64→61→61は保存されたGF結果に限定する。

87次通常多項式模型、monic6関数体模型、cを普通の変数とするblock模型の終了124/14/124を確認した。
元134式保存・凍結入力一致・GB/cofactor成功markerの不在を監査し、unitやeliminantが得られたとはしない。
低次生座標q32=cの追加probeも3GiB上限で失敗し、q32の代数性を証明していない。
全問題の環境未導入とはせず、必要な小さいexact GB/change-matrixは既に成功していることとも区別する。

このpacketから一般有限性、曲面除去、全長5を結論しない。
再開には構造的消去の改善・別の正当化された算法、または明示的な追加計算予算が必要であり、
同じ未変更GBを自動で延長する根拠にはしない。
code/s61_preprocessing_blockers_review_20261103_v62.py と同名 REVIEW report が証拠を保存する。
