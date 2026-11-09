# 現在地 — 導手の可換余法線と平面の全復元

2026-11-09 / v70-conductor-conormal-reconstruction

active queueはresearch/conceptual/GLOBAL_PLAN.md。現行台帳は147成分（rank4:141、rank2:6）・1158採録主張を主張する。
今回その全証明を独立再確認したとはしない。台帳・旧証明書は未変更。

主稿と構造稿は共通TeX四本を読み、既存の全証明を保持する。
Aの標識/modular原始関数と二つの実持上げに加え、Bのcanonical cofactorとPfaffian復元を統合した。
rank三・四の全31場を同じ入口へ接続し、全保存源の全自由係数・残Jacobiをexact入力として保持する。
その全係数恒等式を再生成したが、full150線形kernelの完全性のrank証明は旧入力である。

被約平面五次の全Saito源をA二H+C、B Q+C、C既約Fで共有する。
磁気shearのkernelを因子の留数空間へ同型に写し、ν因子でkernel/cokernelの次元ν−1を証明した。
二因子では非零磁気が23次元軌道、零磁気が22次元境界。既約Fでは全磁気が消え、三因子では商が二次元残る。
六Q+C配置・二H+C・既約Fの全R6とshearの正逆を記号検算した。
個別有限配置・全源帰属・周囲の二次/三次障害は保持し、軌道次元から成分性を推定していない。

新規性未確認。検算はcode/check_quadric_residue_reconstruction.py、入力はcode/data/quadric_cofactor_sources.json。
結果はreports/quadric_residue_checks.json、組版・目視と保存範囲はreports/quadric_residue_validation.json。
Cの導手支持について、rank零・二の可換余法線とflat作用を共通化した。
planeではcanonical liftの全weight源を(45,30,9)、rank零を(35,30,9)として正逆付きで復元した。
restriction・余法線作用を固定してもbinary quadratic三係数が残る。全higherと非斉次Jacobiを保持する。
leading Pfaffianをrestrictionで割ったsectionの平行輸送から、二次判別式と四次不変量を同じ原理で説明する。
rank零可換性、rank二module、既約cubic零因子の全source、m2/m4輸送の重複証明を実際に置換した。
全150係数で制約rank42、Euler24、全源84の完全性と全自由係数恒等式を記号検算した。
新規性未確認。旧有限場合分け・共鳴・monodromy・全支持帰属・成分性は保持し、全147の独立検証ではない。
新検算はcode/check_conductor_conormal_reconstruction.py、結果はreports/conductor_conormal_checks.json。
組版と目視はreports/conductor_conormal_validation.json。次は非plane認識と正規化上の降下・留数を接続する。
GitHub、bundle生成/検証、Library保存は別状態として実結果だけを報告する。
