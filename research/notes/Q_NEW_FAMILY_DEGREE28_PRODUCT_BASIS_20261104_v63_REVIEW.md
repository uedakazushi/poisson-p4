# 全Inv28の積基底：独立係数監査

2026-11-04。producer/nativeを呼ばず、二変数の疎辞書積で全117列・各435係数を再構成した。
入口の80配列は、凍結済みのCD、full Hessian Inv24、ambient Poisson Inv14の各選択列へ
全係数を直接照合した。同じ既存formal planeを使い、新しい点を選んでいない。

F4×Inv24の67本、primitive12×二つのocticの40本、primitive16×primitive12の4本、
Inv14対称積の6本が全てnative列と一致した。選択108×108行列の全11664係数も独立照合し、
mod32003の行列式27379を独自の消去で求めた。
元の16 Cartan重みからcharacterを再計算するとdimInv28=108である。
従ってQQ上の独立108本は全Inv28を張る。標数零の生成性を有限体標本だけで推測していない。

この監査は普遍基底の完成だけを受理する。actual NEW sourceへの28次制限や
定義次数の改善、方程式の構成はまだ行っていない。
