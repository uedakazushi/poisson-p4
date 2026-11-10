# Poisson部分多様体と射影四空間の成分交叉：研究方針

2026-11-10 / v10。再開順はルートRESEARCH_POLICY.md → REPRODUCE_GITHUB.md → 本方針 → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。

## 恒常的原則と保存許可

数学的重要性、新規性、理解の深化、正確性、応用可能性、実現可能性を総合して方針を柔軟に変更する。証明付き／計算検証／未完成／予想／既知・新規性未確認を明示し、必要な仮定、証明の所在、障害を保持する。日本語LuaLaTeXの原稿を正本とし、再取得困難な入力と最小限の証明書を残すが、重複・作業ログ・旧版を蓄積しない。ユーザーの継続保存許可の現行運用では**GitHub mainが唯一の正本**。非force pushとremote先端の再確認が必要で、Library保存は実際のpush失敗時の緊急退避に限る。GitHub Actionsは使わず、PDF・中間生成物は原則追跡しない。既存Libraryファイルを明示指示なく削除しない。

## 現在の成果

- **局所次元と交叉（証明付き、v10）**：q=u0²+2u1w0+2zw1、i_(P0)Ω=dq∧i_ξ(du0∧du1∧dz)、ξ=(1,1,-2)。全ambient Poissonスキームの接空間次元30、局所Krull次元は**正確に28**。二つの全150方向の線形化消滅汎関数を接空間基底へ制限すると互いに素な二次式9*t21*t29と-4*t25*t27。Log(1,1,1,2)源の28次元像は350×28微分の非零28次minorで独立に証明。完全な147成分分類を仮定しない。
- **異なる被約成分との交叉（証明付き）**：P0の対数型二形式はq*u0*u1*z*dlog(u0/u1)∧dlog(q³/(u1⁴z²))。一方、P0+aK1はa≠0で階数四。従ってp=[P0]はrank2対数型の28次元閉包と少なくとも一つのrank4既約成分との交点。rank4成分ID、ambient非被約性は未解決。
- **保持すべき既存結果**：固定Q,σ延長ファイバー Spec C[a,b,c]/(ac,c²)、embedded cの全ambient二次障害、一般smooth Q3のRat/Log/Aff全延長（一点／A¹／二点）、共鳴Log直交性wallと等方的境界三例、smooth cubicの一意延長、一般CI(2,2) del Pezzo剛性、P³横断対称性、既存147成分分類の証明と有限証明書。文献上の全面的優先権・外部査読は未確認。

## 未解決の核心と着手順

最優先は**pを通るrank4既約成分のIDを少なくとも一つ具体的に決定すること**。Pfaffian q*u1²*zを持つPoisson曲線を既存カタログのquadric+hyperplane型・T族の実源と照合し、包含の根拠を退化族で示す。次に全二次接錐の既約分解、固定ファイバーのac関係のambient二次障害、ambient局所環の被約性を研究。単なる局所場合分けになればdel Pezzo例外・P³対称性へ切り替える。

## 原稿と計算

主証明 papers/quadric_log_component_intersection.tex（incidence主稿にinclude）、従来の障害 papers/quadric_ambient_quadratic_obstruction.tex、固定ファイバー papers/quadric_log_isotropic_boundary.tex。新しい厳密producer code/subvariety_incidence/check_local_dimension28.py と保存証明書 research/subvariety_incidence/certificates/local_dimension28.json。手動再現 bash scripts/verify_github_source.sh incidence。生成、計算検証、組版、GitHub保存を別々の状態として報告。
