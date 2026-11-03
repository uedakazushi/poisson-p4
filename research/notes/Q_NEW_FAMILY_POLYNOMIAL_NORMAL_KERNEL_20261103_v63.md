# 二つの普遍linear polynomial法線とgeneric Jacobi rank98の短い証明

2026年11月3日。前の三変数native rank証明を変更せず、独立に短い上界証明を与える。
正しいgeneric二変数normalizationは、既存の全150係数監査済みの**z-h shearを含む**変換を使う。
純粋なdiagonal変換だけでは定数{y,z}が変わるので使わない。
D(n0*n1)上でn1=1へ移る密な族の係数体はQQ(n0,n2)である。

全350×126 physical canonical Jacobi行列を同じ密な族で再構成する。
その多項式kernelをdegree1で探すと、2041×378、9384非零係数の有理線形系から
29個のpolynomial syzygyが得られた。
そのうち二つは既存の群・scalar・二変数実族の26次元像の外にある。
発見方法だけには依存せず、保存する二つの126×2 polynomial vector Tについて、
全350行のJ*T=0を直接確認する。

既存の26個の実際のPoisson tangent方向も合わせて、126×28 polynomial kernel Kを作る。
全350×28積が恒等的に0、(n0,n2)=(1,1)で28×28 minorが非零である。
従ってgeneric kernel dimension≥28、Jacobian rank≤126−28=98。
別に同じ点で原350×126の98×98 minorが非零なので、generic rankは厳密に98となる。
この上界に192秒の有理関数消去の結果や標本rankの小ささは必要ない。
26個の方向が実族に属することは、全Poisson恒等式の微分と実GL/scalar作用による。

正しいsheared PGL normalizationの密性から、これが新しい全成分のgeneric rankでもある。
二つのTはactual26と合わせてgeneric tangent kernel全体を張るので、
generic二法線のobstruction計算に使える。
原物理二次Jacobiの全350×3 Hessian列も保存するが、
そのcokernel rankやgeneric Artin長は本packetだけで結論しない。

* `code/q_new_family_polynomial_normal_kernel_20261103_v63.py`
* `reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63/{report.json,certificate.json.gz}`

全係数・二つのminor・全Hessian・正しいnormalization監査への入力SHA256を保存する。
