# 実Tの一次層と既知二次span：独立Schur・固定rank監査

2026-11-04。実X0/X1、μT=5、complete116×936提示を先行独立peerで固定する。
作者helper・新しいrank探索・source係数・jet・348行列は実行しない。

27個のconstant source basisを消す五つのleft annihilatorを、pivot順の前進消去ではなく
逆順の係数back-substitutionで構成する。各rowが27列を殺し、残る五座標にidentityとなることを
確認して、全741 Q1を再構築する。84 unit wedgeの正規化Q1は、既存の全column factorizationを
三角的に対にして再構築する。全657 residual linear関係をSchur公式と照合する。
保存されたlinear rank証明書の全引算・全column factorization・非零minor積を検査しrank37と確定する。
従って `dim(T/m²T)=5+(40−37)=8`、初め二層HFは(5,3)である。

選択37本の実linear関係の8変数倍を全296列再構築する。
literal d3の一次を、保存済みの84 unit wedge basisとそのbeta座標で再合成して全係数を確認する。
その座標と実normalized Q1から、正しい符号の `−Qstar1 Wstar0^-1 d3_linear` を全168列で検査する。
既知464列は真のdegree2初期関係である。保存されたrank証明書を検査して、この**既知span**のrank177を認める。

180方向を張っていないことは、実m²T/m³Tの次元3や非零性を意味しない。
未入力の他の二次cycle関係が残りを消す可能性がある。
したがってm²T零、全長8、全長>8、K torsion、元B18深さのどれも結論しない。
専用checker `code/ambient_B18_comparison_T_linear_quadratic_review_20261104_v63.py` と同名REVIEWに保存する。
