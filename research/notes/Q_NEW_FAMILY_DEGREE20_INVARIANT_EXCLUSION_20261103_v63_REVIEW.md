# 新first-jet hypersurfaceの20次除外：独立QQ整数word監査

2026-11-03。著者producer/native evaluatorを再実行せず、選択36式・36点の必要十分なexact minorを別engineで確認する。

全36点について、実dense S61 mapの全DAGを独立にFractionでQQ評価し、全inverse node、value pivot、Pfaffianを確認する。
primitive β24と全16 basis係数を再構成し、samples/native-inputへ全scalarをbindingする。
選択点には旧20点の外にある**16個の新しい実QQ source点**が含まれる。
旧点記録を36点へコピーしただけの監査ではない。
generic mapの全parameter恒等式は既存の独立dense-map監査を明示的に使用する。

fresh standard sp4 basisからtrace Gram inverseを取り直し、ω、ω⁻¹、全10表現行列からθ,Hを再構成する。
各点でθ/Hを共通分母で**整数行列**へ変換し、選択36wordを逐次prefix積で計算する。
著者nativeのrational matrix・長さ5+5分割法とは異なる。
全1,296個のselected QQ値が著者証拠と一致し、整数36×36に独立Bareissを施してrow分母を戻す。
正確な有理detは著者の8320/5697-bit determinantに等しく、mod32003は9474で非零。
全3,888個のQQ値を再計算したとは主張しない。

完全な20次上界も別法で確認する。
各**異なるweightとその重複度**をbinomial係数付きで展開し、positive-root denominator四factorを掛けてweight0を取った。
著者の個別weight coin-DPおよび8 Weyl項交代和とは異なる計算法で、Inv20の次元は36。
次数0〜20の全Weyl multiplicityと各Sym^dの総次元binomial(15+d,d)も照合する。
全16 weightのparityから中心−Iはβに−1で作用し、全奇数次数の不変式は零である。

θ/Hの普遍共役性とtrace不変性、実像がSp4-stable既約cone hypersurfaceであり生成式がSp4不変になること、18次以下の除外は既存独立監査を再使用する。
36個の制限独立式と完全上界36から、全20次制限が単射。
従って生成方程式の次数は**22以上**となる。
22次以上のprobe、方程式の存在次数の同定、generic fiber尽くしは行わない。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_degree20_invariant_review_20261103_v63.py --output <new-output-directory>

全fresh QQ点、整数matrix分母、selected値、full exact determinant、character dimensionsを独立review evidenceへ保存する。
120秒・2.2GB、GBなし。
